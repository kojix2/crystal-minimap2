module Minimap2
  # ---------------------------------------------------------------------------
  # Alignment skeleton — port of align.c
  # Applies ksw2 alignment between query and reference using chained anchors
  # as guide points.
  # ---------------------------------------------------------------------------

  # Generate a simple substitution matrix (m×m, match=a, mismatch=b).
  private def self.gen_simple_mat(m : Int32, a : Int32, b : Int32, sc_ambi : Int32, transition : Int32) : Array(Int8)
    if transition != 0 && transition != b
      ksw_gen_ts_mat(m, a, b, transition, sc_ambi)
    else
      ksw_gen_simple_mat(m, a, b, sc_ambi)
    end
  end

  # Append a CIGAR run to a MmExtra, merging with the previous operation if same.
  private def self.append_cigar(ep : MmExtra, cigar : Array(UInt32)) : Nil
    cigar.each do |entry|
      op = entry & 0xf_u32
      len = entry >> 4
      next if len == 0
      push_cigar(ep.cigar, op.to_i32, len.to_i32)
    end
  end

  # Retrieve a segment of the reference sequence (as 4-bit encoded bases).
  private def self.get_ref_seq(mi : MmIdx, rid : Int32, rs : Int32, re : Int32, is_rev : Bool) : Array(UInt8)
    len = re - rs
    buf = Array(UInt8).new(len, 4_u8)
    return buf if len <= 0
    if is_rev
      mi.getseq_rev(rid.to_u32, rs.to_u32, re.to_u32, buf)
    else
      mi.getseq(rid.to_u32, rs.to_u32, re.to_u32, buf)
    end
    buf
  end

  # Encode a query string (ASCII) into 4-bit bases using SEQ_NT4_TABLE.
  def self.encode_seq(s : String | Bytes) : Array(UInt8)
    bytes = s.is_a?(String) ? s.to_slice : s
    result = Array(UInt8).new(bytes.size, 4_u8)
    i = 0
    while i < bytes.size
      result[i] = SEQ_NT4_TABLE[bytes[i].to_i]
      i += 1
    end
    result
  end

  # Build reverse complement of a 4-bit encoded sequence.
  # Delegates to ksw2.seq_rev_comp (shared implementation).
  def self.rev_comp(seq : Array(UInt8) | Slice(UInt8)) : Array(UInt8)
    seq_rev_comp(seq.size, seq)
  end

  # Core alignment between a pair of query and target windows.
  # Fills in r.p with CIGAR and dp scores.
  private def self.align_pair(opt : MmMapOpt, qlen : Int32, qseq : Array(UInt8) | Slice(UInt8),
                              tlen : Int32, tseq : Array(UInt8) | Slice(UInt8),
                              mat : Array(Int8), w : Int32,
                              end_bonus : Int32, zdrop : Int32, flag : Int32) : KswExtz
    if opt.max_sw_mat > 0 && tlen.to_i64 * qlen > opt.max_sw_mat
      ez = KswExtz.new
      ez.zdropped = true
      return ez
    end

    if (opt.flag & F_SPLICE) != 0
      ksw_exts2(qlen, qseq, tlen, tseq, 5, mat,
        opt.q, opt.e, opt.q2, opt.noncan,
        zdrop, end_bonus, opt.junc_bonus, opt.junc_pen, flag)
    elsif opt.q == opt.q2 && opt.e == opt.e2
      ksw_extz2(qlen, qseq, tlen, tseq, 5, mat,
        opt.q, opt.e, w, zdrop, end_bonus, flag)
    else
      ksw_extd2(qlen, qseq, tlen, tseq, 5, mat,
        opt.q, opt.e, opt.q2, opt.e2, w, zdrop, end_bonus, flag)
    end
  end

  # Compute r.mlen, r.blen, r.div from CIGAR (mirrors mm_update_extra).
  private def self.update_extra(r : MmReg1, qseq : Array(UInt8) | Slice(UInt8), tseq : Array(UInt8) | Slice(UInt8),
                                mat : Array(Int8), q : Int32, e : Int32, is_eqx : Bool) : Nil
    ep = r.p
    return unless ep
    r.blen = r.mlen = 0
    r.is_spliced = false
    qoff = 0; toff = 0
    s = 0.0; max_s = 0.0

    ep.cigar.each do |entry|
      op = (entry & 0xf).to_i32
      len = (entry >> 4).to_i32
      case op
      when CIGAR_MATCH
        n_ambi = 0; n_diff = 0
        len.times do |idx|
          cq = qseq[qoff + idx].to_i32
          ct = tseq[toff + idx].to_i32
          if ct > 3 || cq > 3
            n_ambi += 1
          elsif ct != cq
            n_diff += 1
          end
          s += mat[ct * 5 + cq].to_f
          s = 0.0 if s < 0.0
          max_s = s if s > max_s
        end
        r.blen += len - n_ambi
        r.mlen += len - n_ambi - n_diff
        ep.n_ambi += n_ambi.to_u32
        qoff += len; toff += len
      when CIGAR_INS
        len.times { |idx| n_ambi = qseq[qoff + idx] > 3 ? 1 : 0; ep.n_ambi += n_ambi.to_u32; r.blen += 1 - n_ambi }
        s -= q + e
        s = 0.0 if s < 0.0
        qoff += len
      when CIGAR_DEL
        len.times { |idx| n_ambi = tseq[toff + idx] > 3 ? 1 : 0; ep.n_ambi += n_ambi.to_u32; r.blen += 1 - n_ambi }
        s -= q + e
        s = 0.0 if s < 0.0
        toff += len
      when CIGAR_N_SKIP
        r.is_spliced = true
        toff += len
      end
    end

    ep.dp_max = ep.dp_max0 = max_s.to_i32

    # Convert M to =/X if requested
    if is_eqx
      new_cigar = [] of UInt32
      qoff2 = 0; toff2 = 0
      ep.cigar.each do |entry|
        op = (entry & 0xf).to_i32
        len = (entry >> 4).to_i32
        if op == CIGAR_MATCH
          rem = len
          while rem > 0
            l = 0
            while l < rem && qseq[qoff2 + l] == tseq[toff2 + l]
              l += 1
            end
            if l > 0
              new_cigar << (l.to_u32 << 4 | CIGAR_EQ_MATCH.to_u32)
              qoff2 += l; toff2 += l; rem -= l
            end
            l = 0
            while l < rem && qseq[qoff2 + l] != tseq[toff2 + l]
              l += 1
            end
            if l > 0
              new_cigar << (l.to_u32 << 4 | CIGAR_X_MISMATCH.to_u32)
              qoff2 += l; toff2 += l; rem -= l
            end
          end
        else
          new_cigar << entry
          qoff2 += len if op == CIGAR_INS
          toff2 += len if op == CIGAR_DEL || op == CIGAR_N_SKIP
        end
      end
      ep.cigar = new_cigar
    end
  end

  # Align one region r using its seed chain.
  # This is a simplified version of mm_align_skeleton / mm_realign.
  private def self.align_one_reg(opt : MmMapOpt, mi : MmIdx,
                                 qlen : Int32, qseq_fwd : Array(UInt8),
                                 r : MmReg1, a : Array(Mm128)) : Nil
    is_rev = r.rev?
    rid = r.rid
    qs = r.qs; qe = r.qe; rs = r.rs; re = r.re
    qlen_aln = qe - qs; tlen_aln = re - rs
    return if qlen_aln <= 0 || tlen_aln <= 0

    # Retrieve query sequence for this region
    qseq : Array(UInt8) | Slice(UInt8)
    if is_rev
      qseq = rev_comp(Slice.new(qseq_fwd.to_unsafe + qs, qlen_aln))
    else
      qseq = Slice.new(qseq_fwd.to_unsafe + qs, qlen_aln)
    end

    # Retrieve target sequence
    tseq = get_ref_seq(mi, rid, rs, re, false)

    mat = gen_simple_mat(5, opt.a, opt.b, opt.sc_ambi, opt.transition)

    bw = [opt.bw, tlen_aln, qlen_aln].min
    zdrop = opt.zdrop

    ksw_flag = KSW_EZ_EXTZ_ONLY
    ksw_flag |= KSW_EZ_REV_CIGAR
    ksw_flag |= KSW_EZ_RIGHT if (opt.flag & F_SPLICE) != 0

    ez = align_pair(opt, qlen_aln, qseq, tlen_aln, tseq, mat,
      bw, opt.end_bonus, zdrop, ksw_flag)

    return if ez.zdropped? || ez.n_cigar == 0

    # Create MmExtra and fill it
    ep = MmExtra.new
    ep.dp_score = ez.score
    ep.dp_max = ez.score
    ep.dp_max2 = KSW_NEG_INF
    append_cigar(ep, ez.cigar)
    r.p = ep

    is_eqx = (opt.flag & F_EQX) != 0
    update_extra(r, qseq, tseq, mat, opt.q, opt.e, is_eqx)
  end

  # Alignment skeleton: align all regions with CIGAR output.
  # Mirrors mm_align_skeleton().
  def self.align_skeleton(opt : MmMapOpt, mi : MmIdx, qlen : Int32,
                          qstr : String | Array(UInt8),
                          n_regs_ref : Pointer(Int32), regs : Array(MmReg1),
                          a : Array(Mm128)) : Array(MmReg1)
    return regs unless (opt.flag & F_CIGAR) != 0

    qseq_fwd : Array(UInt8)
    if qstr.is_a?(String)
      qseq_fwd = encode_seq(qstr.to_slice)
    else
      qseq_fwd = qstr
    end

    n_regs = n_regs_ref.value
    n_regs.times do |i|
      r = regs[i]
      next if r.cnt <= 0
      align_one_reg(opt, mi, qlen, qseq_fwd, r, a)
    end

    # Filter after alignment
    filter_regs(opt, qlen, n_regs_ref, regs)

    regs
  end

  # Estimate divergence from minimizer matching.
  # Port of mm_est_err() from esterr.c.
  def self.est_err(mi : MmIdx, qlen : Int32, regs : Array(MmReg1),
                   a : Array(Mm128), n_mini_pos : Int32, mini_pos : Array(UInt64)) : Nil
    return if n_mini_pos == 0

    # Compute average k-mer span
    sum_k = 0_u64
    n_mini_pos.times { |i| sum_k += (mini_pos[i] >> 32) & 0xff }
    avg_k = sum_k.to_f / n_mini_pos

    regs.each do |r|
      r.div = -1.0_f32
      next if r.cnt <= 0

      # Get query position for a minimizer, reverting reverse-strand to
      # the forward-strand coordinate.  Mirrors get_for_qpos() in esterr.c.
      get_qpos = ->(idx : Int32) : Int32 do
        entry = a[idx]
        x = u64_to_i32(entry.y)
        q_span = ((entry.y >> 32) & 0xff).to_i32
        if entry.x >> 63 == 1
          qlen - 1 - (x + 1 - q_span)
        else
          x
        end
      end

      # Binary search for the first anchor's query position in mini_pos.
      # For reverse strand, anchors are stored in reverse query order,
      # so the first forward-strand anchor is at a_off + cnt - 1.
      first_idx = r.rev? ? r.a_off + r.cnt - 1 : r.a_off
      first_qpos = get_qpos.call(first_idx)
      st = -1; en = -1
      lo = 0; hi = n_mini_pos - 1
      while lo <= hi
        mid = ((lo.to_u64 + hi) >> 1).to_i32
        y = u64_to_i32(mini_pos[mid])
        if y < first_qpos
          lo = mid + 1
        elsif y > first_qpos
          hi = mid - 1
        else
          st = mid
          break
        end
      end
      if st < 0
        # Could not find matching minimizer — leave div = -1
        next
      end
      en = st

      # Walk through remaining anchors and match them to mini_pos.
      # For reverse strand, iterate anchors backwards (decreasing index).
      l_ref = mi.seq[r.rid].len
      n_match = 1
      k = 1
      j = st + 1
      while j < n_mini_pos && k < r.cnt
        anchor_idx = r.rev? ? r.a_off + r.cnt - 1 - k : r.a_off + k
        x = get_qpos.call(anchor_idx)
        if x == u64_to_i32(mini_pos[j])
          k += 1
          en = j
          n_match += 1
        end
        j += 1
      end

      n_tot = en - st + 1
      # Add boundary checks (same as C)
      n_tot += 1 if r.qs > avg_k.to_i32 && r.rs > avg_k.to_i32
      n_tot += 1 if qlen - r.qs > avg_k.to_i32 && l_ref - r.re > avg_k.to_i32

      r.div = n_match >= n_tot ? 0.0_f32 : (1.0_f32 - (n_match.to_f / n_tot) ** (1.0 / avg_k)).to_f32
    end
  end
end
