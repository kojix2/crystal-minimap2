module Minimap2
  # ---------------------------------------------------------------------------
  # Main mapping pipeline — port of map.c
  # ---------------------------------------------------------------------------

  # X31 string hash (used to hash query names).
  private def self.x31_hash_string(s : String) : UInt32
    h = 0_u32
    s.each_byte { |c| h = h &* 31_u32 &+ c.to_u32 }
    h
  end

  private def self.wang_hash(key : UInt32) : UInt32
    key = (key ^ 61_u32) ^ (key >> 16)
    key = key &+ (key << 3)
    key = key ^ (key >> 4)
    key = key &* 0x27d4eb2d_u32
    key = key ^ (key >> 15)
    key
  end

  # Filter minimizers overlapping with low-complexity (SDUST) regions.
  private def self.dust_minier(mv : Array(Mm128), l_seq : Int32,
                                seq : String, sdust_thres : Int32) : Int32
    return mv.size if sdust_thres <= 0
    dreg = sdust(seq.to_slice, l_seq, sdust_thres, 64)
    u = 0
    k = 0
    mv.each_with_index do |anchor, j|
      q_pos = (anchor.y.to_u32 >> 1).to_i32
      span  = (anchor.x & 0xff).to_i32
      s = q_pos - (span - 1); e = s + span
      while u < dreg.size && dreg[u].to_i32 <= s; u += 1; end
      if u < dreg.size && (dreg[u] >> 32).to_i32 < e
        # Minimizer overlaps with LCR
        l = 0
        v = u
        while v < dreg.size && (dreg[v] >> 32).to_i32 < e
          ss = [s, (dreg[v] >> 32).to_i32].max
          ee = [e, dreg[v].to_i32].min
          l += ee - ss
          v += 1
        end
        mv[k] = mv[j]; k += 1 if l <= span >> 1
      else
        mv[k] = mv[j]; k += 1
      end
    end
    k
  end

  # Collect query minimizers across all segments.
  private def self.collect_minimizers(opt : MmMapOpt, mi : MmIdx,
                                       n_segs : Int32, qlens : Array(Int32),
                                       seqs : Array(String)) : Array(Mm128)
    mv = [] of Mm128
    sum = 0
    n_segs.times do |i|
      n = mv.size
      mm_sketch(seqs[i], qlens[i], mi.w, mi.k, i.to_u32, (mi.flag & I_HPC) != 0, mv)
      # Adjust segment IDs into concatenated query space
      n.upto(mv.size - 1) do |j|
        mv[j] = Mm128.new(mv[j].x, mv[j].y &+ (sum << 1).to_u64)
      end
      if opt.sdust_thres > 0
        new_n = n + dust_minier(mv[n..], qlens[i], seqs[i], opt.sdust_thres)
        mv.delete_at(new_n, mv.size - new_n) if mv.size > new_n
      end
      sum += qlens[i]
    end
    mv
  end

  # Check if a seed should be skipped (diagonal / dual mode).
  private def self.skip_seed?(flag : Int64, r : UInt64, q : MmSeed,
                               qname : String?, qlen : Int32, mi : MmIdx) : {Bool, Bool}
    is_self = false
    if qname && (flag & (F_NO_DIAG | F_NO_DUAL)) != 0
      s = mi.seq[r >> 32]
      cmp = qname <=> s.name
      if (flag & F_NO_DIAG) != 0 && cmp == 0 && s.len.to_i32 == qlen
        return {true, false} if (r.to_u32 >> 1) == (q.q_pos >> 1)
        is_self = (r & 1) == (q.q_pos & 1)
      end
      return {true, false} if (flag & F_NO_DUAL) != 0 && cmp > 0
    end
    if (flag & (F_FOR_ONLY | F_REV_ONLY)) != 0
      if (r & 1) == (q.q_pos & 1)  # forward
        return {true, false} if (flag & F_REV_ONLY) != 0
      else
        return {true, false} if (flag & F_FOR_ONLY) != 0
      end
    end
    {false, is_self}
  end

  # Collect seed hits and build anchor array (mirrors collect_seed_hits).
  private def self.collect_seed_hits(opt : MmMapOpt, max_occ : Int32, mi : MmIdx,
                                      qname : String?, mv : Array(Mm128), qlen : Int32,
                                      use_heap : Bool) : {Array(Mm128), Int32, Array(UInt64)}
    seeds, n_a, rep_len, mini_pos = collect_matches(qlen, max_occ, opt.max_max_occ, opt.occ_dist, mi, mv)

    a = [] of Mm128

    seeds.each do |q|
      q.n.times do |k|
        r = q.cr[k]
        skip, is_self = skip_seed?(opt.flag, r, q, qname, qlen, mi)
        next if skip

        rid = (r >> 32).to_i32
        rpos = (r.to_u32 >> 1).to_i32
        p = Mm128.max

        if (r & 1) == (q.q_pos & 1)  # forward strand
          p = Mm128.new(
            (r & 0xffffffff00000000_u64) | rpos.to_u64,
            q.q_span.to_u64 << 32 | (q.q_pos >> 1).to_u64
          )
        elsif (opt.flag & F_QSTRAND) == 0  # reverse strand, non-qstrand mode
          p = Mm128.new(
            1_u64 << 63 | (r & 0xffffffff00000000_u64) | rpos.to_u64,
            q.q_span.to_u64 << 32 | (qlen - ((q.q_pos >> 1) + 1 - q.q_span.to_i32) - 1).to_u64
          )
        else  # reverse strand, qstrand mode
          s = mi.seq[rid]
          p = Mm128.new(
            1_u64 << 63 | (r & 0xffffffff00000000_u64) | (s.len.to_i32 - (rpos + 1 - q.q_span.to_i32) - 1).to_u64,
            q.q_span.to_u64 << 32 | (q.q_pos >> 1).to_u64
          )
        end

        # Apply segment ID and flags
        p = Mm128.new(p.x, p.y | (q.seg_id.to_u64 << SEED_SEG_SHIFT))
        p = Mm128.new(p.x, p.y | SEED_TANDEM) if q.is_tandem
        p = Mm128.new(p.x, p.y | SEED_SELF) if is_self

        a << p
      end
    end

    radix_sort_128x(a)
    {a, rep_len, mini_pos}
  end

  # Post-chain processing: set parent / secondary / primary flags.
  private def self.chain_post(opt : MmMapOpt, max_chain_gap_ref : Int32, mi : MmIdx,
                               qlen : Int32, n_segs : Int32, qlens : Array(Int32),
                               n_regs_ref : Pointer(Int32), regs : Array(MmReg1),
                               a : Array(Mm128)) : Nil
    return if (opt.flag & F_ALL_CHAINS) != 0
    set_parent(opt.mask_level, opt.mask_len, regs, opt.a * 2 + opt.b,
               (opt.flag & F_HARD_MLEVEL) != 0, opt.alt_drop)
    if n_segs <= 1
      select_sub(opt.pri_ratio, mi.k * 2, opt.best_n, true, (opt.max_gap * 0.8).to_i32,
                 n_regs_ref, regs)
    end
    # multi-segment version not fully implemented here
  end

  # Map a single query fragment (internal core).
  def self.map_frag_core(mi : MmIdx, n_segs : Int32, qlens : Array(Int32),
                          seqs : Array(String), n_regs_arr : Array(Int32),
                          regs_arr : Array(Array(MmReg1)),
                          opt : MmMapOpt, qname : String? = nil) : Nil
    qlen_sum = qlens.sum
    n_segs.times { |i| n_regs_arr[i] = 0; regs_arr[i] = [] of MmReg1 }
    return if qlen_sum == 0 || n_segs <= 0

    # Compute query hash
    hash = qname && (opt.flag & F_NO_HASH_NAME) == 0 ? x31_hash_string(qname) : 0_u32
    hash ^= wang_hash(qlen_sum.to_u32) &+ wang_hash(opt.seed.to_u32)
    hash = wang_hash(hash)

    # Collect minimizers
    mv = collect_minimizers(opt, mi, n_segs, qlens, seqs)
    seed_mz_flt(mv, opt.mid_occ, opt.q_occ_frac) if opt.q_occ_frac > 0.0_f32

    use_heap = (opt.flag & F_HEAP_SORT) != 0
    a, rep_len, mini_pos = collect_seed_hits(opt, opt.mid_occ, mi, qname, mv, qlen_sum, use_heap)

    is_splice   = (opt.flag & F_SPLICE) != 0
    is_sr       = (opt.flag & F_SR) != 0

    max_chain_gap_qry = is_sr ? [qlen_sum, opt.max_gap].max : opt.max_gap
    max_chain_gap_ref = if opt.max_gap_ref > 0
      opt.max_gap_ref
    elsif opt.max_frag_len > 0
      [opt.max_frag_len - qlen_sum, opt.max_gap].max
    else
      opt.max_gap
    end

    chn_pen_gap  = opt.chain_gap_scale * 0.01_f32 * mi.k
    chn_pen_skip = opt.chain_skip_scale * 0.01_f32 * mi.k

    u = [] of UInt64
    a = if (opt.flag & F_RMQ) != 0
      lchain_rmq(opt.max_gap, opt.rmq_inner_dist, opt.bw, opt.max_chain_skip,
                 opt.rmq_size_cap, opt.min_cnt, opt.min_chain_score,
                 chn_pen_gap, chn_pen_skip, a, u)
    else
      lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt.bw, opt.max_chain_skip,
                opt.max_chain_iter, opt.min_cnt, opt.min_chain_score,
                chn_pen_gap, chn_pen_skip, is_splice, n_segs, a, u)
    end

    # Re-chain with long bw if needed
    n_regs0 = u.size
    if opt.bw_long > opt.bw && (opt.flag & (F_SPLICE | F_SR | F_NO_LJOIN)) == 0 &&
       n_segs == 1 && n_regs0 > 1 && !a.empty?
      st = u64_to_i32(a[0].y); en = u64_to_i32(a[u64_to_i32(u[0] & 0xffffffff_u64) - 1].y)
      if qlen_sum - (en - st) > opt.rmq_rescue_size || en - st > qlen_sum * opt.rmq_rescue_ratio
        n_a_new = u.sum { |uu| u64_to_i32(uu) }
        u.clear
        a = a.first(n_a_new)
        radix_sort_128x(a)
        a = lchain_rmq(opt.max_gap, opt.rmq_inner_dist, opt.bw_long, opt.max_chain_skip,
                       opt.rmq_size_cap, opt.min_cnt, opt.min_chain_score,
                       chn_pen_gap, chn_pen_skip, a, u)
        n_regs0 = u.size
      end
    elsif opt.max_occ > opt.mid_occ && rep_len > 0 && (opt.flag & F_RMQ) == 0
      # Re-chain with higher max_occ for short reads
      rechain = n_regs0 == 0
      if !rechain && n_regs0 > 0
        # Check if best chain covers all segments
        max_sc  = 0; max_i = 0; max_off = 0; off = 0
        n_regs0.times do |i|
          sc = (u[i] >> 32).to_i32
          if sc > max_sc; max_sc = sc; max_i = i; max_off = off; end
          off += u[i].to_i32
        end
        n_chained = 1
        1.upto(u[max_i].to_i32 - 1) do |j|
          if (a[max_off + j].y & SEED_SEG_MASK) != (a[max_off + j - 1].y & SEED_SEG_MASK)
            n_chained += 1
          end
        end
        rechain = n_chained < n_segs
      end
      if rechain
        a, rep_len, mini_pos = collect_seed_hits(opt, opt.max_occ, mi, qname, mv, qlen_sum, use_heap)
        u.clear
        a = lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt.bw, opt.max_chain_skip,
                      opt.max_chain_iter, opt.min_cnt, opt.min_chain_score,
                      chn_pen_gap, chn_pen_skip, is_splice, n_segs, a, u)
        n_regs0 = u.size
      end
    end

    regs0 = gen_regs(hash, qlen_sum, u, a, (opt.flag & F_QSTRAND) != 0)
    if mi.n_alt > 0
      mark_alt(mi, regs0)
      n_regs0_tmp = n_regs0
      hit_sort(pointerof(n_regs0_tmp), regs0, opt.alt_drop)
      n_regs0 = n_regs0_tmp
      regs0 = regs0.first(n_regs0)
    end

    chain_post(opt, max_chain_gap_ref, mi, qlen_sum, n_segs, qlens,
               pointerof(n_regs0), regs0, a)

    unless is_sr || (opt.flag & F_QSTRAND) != 0
      est_err(mi, qlen_sum, regs0, a, mini_pos.size, mini_pos)
      n_regs0 = filter_strand_retained(regs0)
    end

    if n_segs == 1
      n_ref = n_regs0
      regs0 = align_skeleton(opt, mi, qlens[0], seqs[0], pointerof(n_ref), regs0, a)
      if (opt.flag & F_ALL_CHAINS) == 0
        set_parent(opt.mask_level, opt.mask_len, regs0, opt.a * 2 + opt.b,
                   (opt.flag & F_HARD_MLEVEL) != 0, opt.alt_drop)
        n_ref_ptr = n_ref
        select_sub(opt.pri_ratio, mi.k * 2, opt.best_n, false, 0,
                   pointerof(n_ref_ptr), regs0)
        n_ref = n_ref_ptr
        set_sam_pri(regs0)
      end
      set_mapq2(regs0, n_ref, opt.min_chain_score, opt.a, rep_len, is_sr, is_splice)
      n_regs_arr[0] = n_ref
      regs_arr[0] = regs0.first(n_ref)
    else
      segs = seg_gen(hash, n_segs, qlens, regs0, a)
      n_segs.times do |s|
        seg_regs = gen_regs(hash, qlens[s], segs[s].u, segs[s].a, false)
        seg_regs.each { |r| r.seg_split = true; r.seg_id = s.to_u32 }
        n_ref_seg = seg_regs.size
        seg_regs = align_skeleton(opt, mi, qlens[s], seqs[s],
                                  pointerof(n_ref_seg), seg_regs, segs[s].a)
        set_mapq2(seg_regs, n_ref_seg, opt.min_chain_score, opt.a, rep_len, is_sr, is_splice)
        n_regs_arr[s] = n_ref_seg
        regs_arr[s] = seg_regs.first(n_ref_seg)
      end
    end
  end

  # Map a single query string against the index.
  # Mirrors mm_map().
  def self.map(mi : MmIdx, qlen : Int32, seq : String,
               opt : MmMapOpt, qname : String? = nil) : Array(MmReg1)
    n_regs = [0]
    regs   = [[] of MmReg1]
    map_frag_core(mi, 1, [qlen], [seq], n_regs, regs, opt, qname)
    regs[0].first(n_regs[0])
  end

  # Map multiple query sequences (e.g. for fragment mode).
  def self.map_frag(mi : MmIdx, seqs : Array(String), qlens : Array(Int32),
                    opt : MmMapOpt, qname : String? = nil) : Array(Array(MmReg1))
    n_segs  = seqs.size
    n_regs  = Array(Int32).new(n_segs, 0)
    regs    = Array(Array(MmReg1)).new(n_segs) { [] of MmReg1 }
    map_frag_core(mi, n_segs, qlens, seqs, n_regs, regs, opt, qname)
    n_segs.times.map { |i| regs[i].first(n_regs[i]) }.to_a
  end

  # ---------------------------------------------------------------------------
  # High-level Aligner convenience class
  # ---------------------------------------------------------------------------
  class Aligner
    getter idx_opt  : MmIdxOpt
    getter map_opt  : MmMapOpt
    getter indices  : Array(MmIdx)

    def initialize(@idx_opt : MmIdxOpt, @map_opt : MmMapOpt, @indices : Array(MmIdx))
    end

    # Build an Aligner from sequence strings.
    def self.from_strings(seqs : Array(String), names : Array(String)? = nil,
                          preset : String? = nil) : Aligner
      iopt = MmIdxOpt.new; mopt = MmMapOpt.new
      Minimap2.set_opt(nil, iopt, mopt)
      Minimap2.set_opt(preset, iopt, mopt) if preset
      idx = MmIdx.from_strings(iopt.w, iopt.k, (iopt.flag & I_HPC) != 0,
                               iopt.bucket_bits, seqs, names)
      Minimap2.mapopt_update(mopt, idx)
      new(iopt, mopt, [idx])
    end

    # Build an Aligner from a FASTA/FASTQ file.
    def self.from_file(fn : String, preset : String? = nil,
                       n_threads : Int32 = 1) : Aligner?
      iopt = MmIdxOpt.new; mopt = MmMapOpt.new
      Minimap2.set_opt(nil, iopt, mopt)
      Minimap2.set_opt(preset, iopt, mopt) if preset
      reader = MmIdxReader.new(fn, iopt)
      indices = [] of MmIdx
      while (mi = reader.read(n_threads))
        Minimap2.mapopt_update(mopt, mi)
        indices << mi
      end
      reader.close
      return nil if indices.empty?
      new(iopt, mopt, indices)
    end

    # Map a query string; yields MmReg1 hits.
    def map(seq : String, qname : String? = nil, & : MmReg1 -> Nil) : Nil
      @indices.each do |mi|
        regs = Minimap2.map(mi, seq.size, seq, @map_opt, qname)
        regs.each { |r| yield r }
      end
    end

    # Map a query string; returns Array of hits.
    def map(seq : String, qname : String? = nil) : Array(MmReg1)
      result = [] of MmReg1
      map(seq, qname) { |r| result << r }
      result
    end

    # Return the index sequences.
    def seq_names : Array(String)
      @indices.flat_map { |mi| mi.seq.map(&.name) }
    end

    # Return the index sequence lengths.
    def seq_lengths : Array(UInt32)
      @indices.flat_map { |mi| mi.seq.map(&.len) }
    end
  end
end
