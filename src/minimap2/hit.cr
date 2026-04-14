module Minimap2
  # ---------------------------------------------------------------------------
  # Hit processing — port of hit.c
  # ---------------------------------------------------------------------------

  # Simple 64-bit hash (used to compute hash of chain anchor for randomisation)
  private def self.hash64_hit(key : UInt64) : UInt64
    key = (~key &+ (key << 21))
    key = key ^ (key >> 24)
    key = (key &+ (key << 3)) &+ (key << 8)
    key = key ^ (key >> 14)
    key = (key &+ (key << 2)) &+ (key << 4)
    key = key ^ (key >> 28)
    key = key &+ (key << 31)
    key
  end

  # Compute fuzzy mlen/blen from chain anchors.
  private def self.cal_fuzzy_len(r : MmReg1, a : Array(Mm128)) : Nil
    r.mlen = r.blen = 0
    return if r.cnt <= 0
    r.mlen = r.blen = ((a[r.a_off].y >> 32) & 0xff).to_i32
    (r.a_off + 1).upto(r.a_off + r.cnt - 1) do |i|
      span = ((a[i].y >> 32) & 0xff).to_i32
      tl = u64_to_i32(a[i].x) - u64_to_i32(a[i - 1].x)
      ql = u64_to_i32(a[i].y) - u64_to_i32(a[i - 1].y)
      r.blen += tl > ql ? tl : ql
      r.mlen += if tl > span && ql > span
                  span
                elsif tl < ql
                  tl
                else
                  ql
                end
    end
  end

  # Fill coordinates of a region from its chain anchors.
  private def self.reg_set_coor(r : MmReg1, qlen : Int32, a : Array(Mm128), is_qstrand : Bool) : Nil
    k = r.a_off
    q_span = ((a[k].y >> 32) & 0xff).to_i32
    r.rev = (a[k].x >> 63) == 1
    r.rid = ((a[k].x << 1) >> 33).to_i32
    rs_cand = u64_to_i32(a[k].x) + 1 - q_span
    r.rs = rs_cand > 0 ? rs_cand : 0
    r.re = u64_to_i32(a[k + r.cnt - 1].x) + 1
    if !r.rev? || is_qstrand
      r.qs = u64_to_i32(a[k].y) + 1 - q_span
      r.qe = u64_to_i32(a[k + r.cnt - 1].y) + 1
    else
      r.qs = qlen - (u64_to_i32(a[k + r.cnt - 1].y) + 1)
      r.qe = qlen - (u64_to_i32(a[k].y) + 1 - q_span)
    end
    cal_fuzzy_len(r, a)
  end

  # Convert chains (u[], a[]) to an Array(MmReg1), sorted by descending score.
  # Mirrors mm_gen_regs().
  def self.gen_regs(hash : UInt32, qlen : Int32, u : Array(UInt64), a : Array(Mm128), is_qstrand : Bool) : Array(MmReg1)
    n_u = u.size
    return [] of MmReg1 if n_u == 0

    # Build sort keys: (score ^ hash) | anchor_offset
    z = Array(Mm128).new(n_u, Mm128.max)
    k = 0
    n_u.times do |i|
      h = (hash64_hit(hash64_hit(a[k].x) &+ hash64_hit(a[k].y) ^ hash.to_u64) & 0xffffffff_u64).to_u32
      z[i] = Mm128.new((u[i] ^ h.to_u64), (k.to_u64 << 32) | (u[i] & 0xffffffff_u64))
      k += u64_to_i32(u[i])
    end
    radix_sort_128x(z)

    # Reverse so largest score comes first
    z.reverse!

    regs = Array(MmReg1).new(n_u) { MmReg1.new }
    n_u.times do |i|
      ri = regs[i]
      ri.id = i
      ri.parent = PARENT_UNSET
      ri.score = ri.score0 = u64_to_i32(z[i].x >> 32)
      ri.hash = (z[i].x & 0xffffffff_u64).to_u32
      ri.cnt = u64_to_i32(z[i].y)
      ri.a_off = (z[i].y >> 32).to_i32
      ri.div = -1.0_f32
      reg_set_coor(ri, qlen, a, is_qstrand)
    end
    regs
  end

  # Mark ALT sequences in hits.
  def self.mark_alt(mi : MmIdx, regs : Array(MmReg1)) : Nil
    return if mi.n_alt == 0
    regs.each do |r|
      r.is_alt = (mi.seq[r.rid].is_alt != 0)
    end
  end

  private def self.alt_score(score : Int32, alt_diff_frac : Float32) : Int32
    return score if score < 0
    s = (score * (1.0 - alt_diff_frac) + 0.499).to_i32
    s > 0 ? s : 1
  end

  # Sort hits by score (descending), removing cnt==0 entries.
  # Mirrors mm_hit_sort().
  def self.hit_sort(n_regs_ref : Pointer(Int32), regs : Array(MmReg1), alt_diff_frac : Float32) : Nil
    n = n_regs_ref.value
    return if n <= 1

    # Keep only hits with cnt > 0 or inv flag set
    aux = Array(Mm128).new
    valid = Array(MmReg1).new

    n.times do |i|
      r = regs[i]
      if r.inv? || r.cnt > 0
        score = (ep = r.p) ? ep.dp_max : r.score
        score = alt_score(score, alt_diff_frac) if r.is_alt?
        aux << Mm128.new((score.to_u64 << 32) | r.hash.to_u64, i.to_u64)
      end
    end

    radix_sort_128x(aux)
    n_aux = aux.size
    n_aux.times do |i|
      valid << regs[(aux[n_aux - 1 - i].y).to_i32]
    end
    n_aux.times { |i| regs[i] = valid[i] }
    n_regs_ref.value = n_aux
  end

  # Set SAM primary flag; returns number of primary hits.
  def self.mark_sam_pri(regs : Array(MmReg1)) : Int32
    n_pri = 0
    regs.each do |r|
      if r.id == r.parent
        n_pri += 1
        r.sam_pri = (n_pri == 1)
      else
        r.sam_pri = false
      end
    end
    n_pri
  end

  # Keep id/parent in sync after hits are removed.
  def self.sync_regs(regs : Array(MmReg1)) : Nil
    return if regs.empty?
    max_id = regs.max_of(&.id)
    return if max_id < 0
    tmp = Array(Int32).new(max_id + 1, -1)
    regs.each_with_index { |r, i| tmp[r.id] = i if r.id >= 0 }
    regs.each_with_index do |r, i|
      r.id = i
      r.parent = if r.parent == PARENT_TMP_PRI
                   i
                 elsif r.parent >= 0 && tmp[r.parent] >= 0
                   tmp[r.parent]
                 else
                   PARENT_UNSET
                 end
    end
    mark_sam_pri(regs)
  end

  # Set parent/secondary relationships.
  # Mirrors mm_set_parent().
  def self.set_parent(mask_level : Float32, mask_len : Int32,
                      regs : Array(MmReg1), sub_diff : Int32,
                      hard_mask_level : Bool, alt_diff_frac : Float32) : Nil
    n = regs.size
    return if n == 0
    n.times { |i| regs[i].id = i }
    w = [] of Int32 # primary hit indices
    w << 0
    regs[0].parent = 0

    (1...n).each do |i|
      ri = regs[i]
      si = ri.qs; ei = ri.qe
      uncov = 0
      cov = [] of UInt64

      unless hard_mask_level
        w.each do |wj|
          rp = regs[wj]
          sj = rp.qs; ej = rp.qe
          next if ej <= si || sj >= ei
          sj = si if sj < si
          ej = ei if ej > ei
          cov << (sj.to_u64 << 32 | ej.to_u64)
        end
        if cov.size > 0
          radix_sort_64(cov)
          x = si
          cov.each do |cv|
            uncov += (cv >> 32).to_i32 - x if (cv >> 32).to_i32 > x
            x = [u64_to_i32(cv), x].max
          end
          uncov += ei - x if ei > x
        end
      end

      # Find a parent
      parent_found = false
      w.each do |wj|
        rp = regs[wj]
        sj = rp.qs; ej = rp.qe
        next if ej <= si || sj >= ei
        min = [ej - sj, ei - si].min
        max = [ej - sj, ei - si].max
        ol = [0, [si < sj ? ei - sj : ej - si, min].min].max
        if ol.to_f / min - uncov.to_f / max > mask_level && uncov <= mask_len
          ri.parent = rp.parent
          sci = ri.score
          sci = alt_score(sci, alt_diff_frac) if !rp.is_alt? && ri.is_alt?
          rp.subsc = [rp.subsc, sci].max
          cnt_sub = ri.cnt >= rp.cnt ? 1 : 0
          if (ep = rp.p) && (ei2 = ri.p)
            if rp.rid != ri.rid || rp.rs != ri.rs || rp.re != ri.re || ol != min
              sci2 = ei2.dp_max
              sci2 = alt_score(sci2, alt_diff_frac) if !rp.is_alt? && ri.is_alt?
              ep.dp_max2 = [ep.dp_max2, sci2].max
              cnt_sub = 1 if ep.dp_max - ei2.dp_max <= sub_diff
            end
          end
          rp.n_sub += cnt_sub
          parent_found = true
          break
        end
      end

      unless parent_found
        w << i
        ri.parent = i
        ri.n_sub = 0
      end
    end
  end

  # Select secondary hits up to best_n.
  def self.select_sub(pri_ratio : Float32, min_diff : Int32, best_n : Int32,
                      check_strand : Bool, min_strand_sc : Int32,
                      n_regs_ref : Pointer(Int32), regs : Array(MmReg1)) : Nil
    return if pri_ratio <= 0.0_f32 || n_regs_ref.value == 0
    n = n_regs_ref.value
    k = 0
    n_2nd = 0

    n.times do |i|
      r = regs[i]
      p = r.parent
      if p == i || r.inv?
        regs[k] = r; k += 1
      elsif (r.score >= regs[p].score * pri_ratio || r.score + min_diff >= regs[p].score) && n_2nd < best_n
        rp = regs[p]
        unless r.qs == rp.qs && r.qe == rp.qe && r.rid == rp.rid && r.rs == rp.rs && r.re == rp.re
          regs[k] = r; k += 1; n_2nd += 1
        end
      elsif check_strand && n_2nd < best_n && r.score > min_strand_sc && r.rev? != regs[p].rev?
        r.strand_retained = true
        regs[k] = r; k += 1; n_2nd += 1
      end
    end

    if k != n
      regs.delete_at(k, regs.size - k) if regs.size > k # remove excess
      trimmed = regs.first(k)
      regs.replace(trimmed)
      sync_regs(regs)
    end
    n_regs_ref.value = k
  end

  # Filter out low-quality hits.
  def self.filter_regs(opt : MmMapOpt, qlen : Int32,
                       n_regs_ref : Pointer(Int32), regs : Array(MmReg1)) : Nil
    k = 0
    n_regs_ref.value.times do |i|
      r = regs[i]
      flt = false
      if !r.inv? && !r.seg_split? && r.cnt < opt.min_cnt
        flt = true
      end
      if ep = r.p
        flt = true if r.mlen < opt.min_chain_score
        flt = true if ep.dp_max < opt.min_dp_max
        flt = true if r.qs > qlen * opt.max_clip_ratio && qlen - r.qe > qlen * opt.max_clip_ratio
      end
      unless flt
        regs[k] = r; k += 1
      end
    end
    n_regs_ref.value = k
  end

  # Filter hits where strand_retained quality is poor.
  def self.filter_strand_retained(regs : Array(MmReg1)) : Int32
    k = 0
    regs.each_with_index do |r, i|
      p = r.parent
      if !r.strand_retained? || p < 0 || p >= regs.size ||
         r.div < regs[p].div * 5.0_f32 || r.div < 0.01_f32
        regs[k] = regs[i]; k += 1
      end
    end
    k
  end

  # Compact the anchor array to only include anchors referenced by regs.
  def self.squeeze_a(regs : Array(MmReg1), a : Array(Mm128)) : Int32
    n = regs.size
    aux = Array(UInt64).new(n) { |i| (regs[i].a_off.to_u64 << 32) | i.to_u64 }
    radix_sort_64(aux)
    as_pos = 0
    aux.each do |av|
      r = regs[(av & 0xffffffff_u64).to_i32]
      if r.a_off != as_pos
        r.cnt.times { |j| a[as_pos + j] = a[r.a_off + j] }
        r.a_off = as_pos
      end
      as_pos += r.cnt
    end
    as_pos
  end

  # Split region r at anchor position n, creating r2 as the second half.
  def self.split_reg(r : MmReg1, r2 : MmReg1, n : Int32, qlen : Int32,
                     a : Array(Mm128), is_qstrand : Bool) : Nil
    return if n <= 0 || n >= r.cnt
    r2.id = -1
    r2.sam_pri = false
    r2.p = nil
    r2.split_inv = false
    r2.cnt = r.cnt - n
    r2.score = (r.score * (r2.cnt.to_f / r.cnt) + 0.499).to_i32
    r2.a_off = r.a_off + n
    r2.parent = r.parent == r.id ? PARENT_TMP_PRI : r.parent
    reg_set_coor(r2, qlen, a, is_qstrand)
    r.cnt -= r2.cnt
    r.score -= r2.score
    reg_set_coor(r, qlen, a, is_qstrand)
    r.split |= 1_u32; r2.split |= 2_u32
  end

  # Set MAPQ scores for hits.
  # Simplified version of mm_set_mapq2 from hit.c.
  def self.set_mapq2(regs : Array(MmReg1), n_regs : Int32,
                     min_chain_sc : Int32, match_sc : Int32,
                     rep_len : Int32, is_sr : Bool, is_splice : Bool) : Nil
    n_regs.times do |i|
      r = regs[i]
      next if r.cnt == 0
      if ep = r.p
        # Estimate mapq from dp_max and dp_max2
        if r.parent == i # primary
          if ep.dp_max2 < 0
            r.mapq = 60_u32
          else
            mapq = (40.0 * (1.0 - ep.dp_max2.to_f / ep.dp_max) * Math.log(r.cnt.to_f + 1) + 0.499).to_u32
            mapq = [mapq, 60_u32].min
            r.mapq = mapq
          end
        else
          r.mapq = 0_u32
        end
      else
        # No alignment, use chaining score
        if r.parent == i
          r.mapq = 60_u32
        else
          r.mapq = 0_u32
        end
      end
    end
  end

  # Generate per-segment chains from the full (concatenated) chain set.
  # Mirrors mm_seg_gen().
  def self.seg_gen(hash : UInt32, n_segs : Int32, qlens : Array(Int32),
                   regs0 : Array(MmReg1), a : Array(Mm128)) : Array(MmSeg)
    acc_qlen = Array(Int32).new(n_segs + 1, 0)
    (1..n_segs).each { |s| acc_qlen[s] = acc_qlen[s - 1] + qlens[s - 1] }
    qlen_sum = acc_qlen[n_segs]

    segs = Array(MmSeg).new(n_segs) { MmSeg.new }
    n_regs0 = regs0.size

    # Initialize per-segment u arrays with scores
    n_segs.times do |s|
      segs[s].u = Array(UInt64).new(n_regs0) { |i| regs0[i].score.to_u64 << 32 }
    end

    # Count anchors per segment per chain
    n_regs0.times do |i|
      r = regs0[i]
      r.cnt.times do |j|
        sid = ((a[r.a_off + j].y & SEED_SEG_MASK) >> SEED_SEG_SHIFT).to_i32
        segs[sid].u[i] += 1
        segs[sid].n_a += 1
      end
    end

    # Compact u arrays
    n_segs.times do |s|
      sr = segs[s]
      new_u = [] of UInt64
      sr.u.each { |v| new_u << v if v.to_i32 != 0 }
      sr.u = new_u
      sr.n_u = new_u.size
      sr.a = Array(Mm128).new(sr.n_a, Mm128.max)
      sr.n_a = 0
    end

    # Fill per-segment anchors
    n_regs0.times do |i|
      r = regs0[i]
      r.cnt.times do |j|
        a1 = a[r.a_off + j]
        sid = ((a1.y & SEED_SEG_MASK) >> SEED_SEG_SHIFT).to_i32
        shift = (a1.x >> 63) != 0 ? (qlen_sum - (qlens[sid] + acc_qlen[sid])).to_u64 : acc_qlen[sid].to_u64
        a1 = Mm128.new(a1.x, a1.y &- shift)
        segs[sid].a[segs[sid].n_a] = a1
        segs[sid].n_a += 1
      end
    end

    segs
  end
end
