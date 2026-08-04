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

  private def self.seq_reverse(seq : Array(UInt8)) : Nil
    i = 0
    j = seq.size - 1
    while i < j
      seq[i], seq[j] = seq[j], seq[i]
      i += 1
      j -= 1
    end
  end

  private def self.seq_window(seq : Array(UInt8), st : Int32, en : Int32) : Array(UInt8)
    len = en - st
    return [] of UInt8 if len <= 0
    Array(UInt8).new(len) { |i| seq[st + i] }
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

  private def self.anchor_rid(a : Mm128) : Int32
    ((a.x << 1) >> 33).to_i32
  end

  private def self.anchor_rev?(a : Mm128) : Bool
    (a.x >> 63) != 0
  end

  private def self.anchor_span(a : Mm128) : Int32
    ((a.y >> 32) & 0xff).to_i32
  end

  private def self.adjust_minier(mi : MmIdx, qseq0 : Array(Array(UInt8)), a : Mm128) : {Int32, Int32}
    if (mi.flag & I_HPC) != 0
      qseq = qseq0[anchor_rev?(a) ? 1 : 0]
      q = u64_to_i32(a.y)
      c = qseq[q]
      i = q - 1
      while i > 0 && qseq[i] == c
        i -= 1
      end
      r = u64_to_i32(a.x) + 1 - get_hplen_back(mi, anchor_rid(a), u64_to_i32(a.x))
      {r, i + 1}
    else
      {u64_to_i32(a.x) - (mi.k >> 1), u64_to_i32(a.y) - (mi.k >> 1)}
    end
  end

  private def self.get_hplen_back(mi : MmIdx, rid : Int32, x : Int32) : Int32
    off0 = mi.seq[rid].offset.to_i64
    off = off0 + x
    c = mi.seq4_get(off.to_u64)
    i = off - 1
    while i >= off0
      break if mi.seq4_get(i.to_u64) != c
      i -= 1
    end
    (off - i).to_i32
  end

  private def self.fix_cigar(r : MmReg1, qseq : Array(UInt8) | Slice(UInt8),
                             tseq : Array(UInt8) | Slice(UInt8)) : {Int32, Int32}
    ep = r.p
    return {0, 0} unless ep
    return {0, 0} if ep.cigar.size <= 1

    toff = 0
    qoff = 0
    to_shrink = false
    k = 0
    while k < ep.cigar.size
      op = (ep.cigar[k] & 0xf).to_i32
      len = (ep.cigar[k] >> 4).to_i32
      to_shrink = true if len == 0
      case op
      when CIGAR_MATCH
        toff += len
        qoff += len
      when CIGAR_INS, CIGAR_DEL
        if k > 0 && k < ep.cigar.size - 1 &&
           (ep.cigar[k - 1] & 0xf) == CIGAR_MATCH &&
           (ep.cigar[k + 1] & 0xf) == CIGAR_MATCH
          prev_len = (ep.cigar[k - 1] >> 4).to_i32
          l = 0
          if op == CIGAR_INS
            while l < prev_len && qoff - 1 - l >= 0 && qoff + len - 1 - l < qseq.size &&
                  qseq[qoff - 1 - l] == qseq[qoff + len - 1 - l]
              l += 1
            end
          else
            while l < prev_len && toff - 1 - l >= 0 && toff + len - 1 - l < tseq.size &&
                  tseq[toff - 1 - l] == tseq[toff + len - 1 - l]
              l += 1
            end
          end
          if l > 0
            ep.cigar[k - 1] = ep.cigar[k - 1] &- (l.to_u32 << 4)
            ep.cigar[k + 1] = ep.cigar[k + 1] &+ (l.to_u32 << 4)
            qoff -= l
            toff -= l
          end
          to_shrink = true if l == prev_len
        end
        if op == CIGAR_INS
          qoff += len
        else
          toff += len
        end
      when CIGAR_N_SKIP
        toff += len
      end
      k += 1
    end

    k = 0
    while k < ep.cigar.size - 2
      op0 = (ep.cigar[k] & 0xf).to_i32
      op1 = (ep.cigar[k + 1] & 0xf).to_i32
      if op0 > 0 && op0 + op1 == 3
        s_ins = 0_u32
        s_del = 0_u32
        l = k
        while l < ep.cigar.size
          op = (ep.cigar[l] & 0xf).to_i32
          len = ep.cigar[l] >> 4
          if op == CIGAR_INS || op == CIGAR_DEL || len == 0
            s_ins += len if op == CIGAR_INS
            s_del += len if op == CIGAR_DEL
          else
            break
          end
          l += 1
        end
        if s_ins > 0 && s_del > 0 && l - k > 2
          ep.cigar[k] = (s_ins << 4) | CIGAR_INS.to_u32
          ep.cigar[k + 1] = (s_del << 4) | CIGAR_DEL.to_u32
          m = k + 2
          while m < l
            ep.cigar[m] &= 0xf_u32
            m += 1
          end
          to_shrink = true
        end
        k = l
      else
        k += 1
      end
    end

    if to_shrink
      squeezed = [] of UInt32
      ep.cigar.each do |entry|
        next if (entry >> 4) == 0
        op = entry & 0xf
        if !squeezed.empty? && (squeezed[-1] & 0xf) == op
          squeezed[-1] = squeezed[-1] &+ ((entry >> 4) << 4)
        else
          squeezed << entry
        end
      end
      ep.cigar = squeezed
    end

    qshift = 0
    tshift = 0
    if !ep.cigar.empty?
      op = (ep.cigar[0] & 0xf).to_i32
      len = (ep.cigar[0] >> 4).to_i32
      if op == CIGAR_INS || op == CIGAR_DEL
        if op == CIGAR_INS
          if r.rev?
            r.qe -= len
          else
            r.qs += len
          end
          qshift = len
        else
          r.rs += len
          tshift = len
        end
        ep.cigar.delete_at(0)
      end
    end
    {qshift, tshift}
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
    ep.n_ambi = 0_u32
    r.blen = r.mlen = 0
    r.is_spliced = false
    qoff = 0; toff = 0
    s = 0.0; max_s = 0.0

    qshift, tshift = fix_cigar(r, qseq, tseq)
    qoff += qshift
    toff += tshift

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
            while l < rem && qseq[qshift + qoff2 + l] == tseq[tshift + toff2 + l]
              l += 1
            end
            if l > 0
              new_cigar << (l.to_u32 << 4 | CIGAR_EQ_MATCH.to_u32)
              qoff2 += l; toff2 += l; rem -= l
            end
            l = 0
            while l < rem && qseq[qshift + qoff2 + l] != tseq[tshift + toff2 + l]
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
    if (opt.flag & (F_SPLICE | F_SR | F_QSTRAND)) == 0 && (mi.flag & I_HPC) == 0 && r.cnt > 0
      align_one_reg_guided(opt, mi, qlen, qseq_fwd, r, a)
      return if r.p
    end

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

  private def self.align_one_reg_guided(opt : MmMapOpt, mi : MmIdx,
                                        qlen : Int32, qseq_fwd : Array(UInt8),
                                        r : MmReg1, a : Array(Mm128)) : Nil
    return if r.cnt <= 0

    qseq_rev = rev_comp(qseq_fwd)
    qseq0 = [qseq_fwd, qseq_rev]
    rev = r.rev?
    qseq_all = qseq0[rev ? 1 : 0]
    rid = r.rid
    ctg_len = mi.seq[rid].len.to_i32
    mat = gen_simple_mat(5, opt.a, opt.b, opt.sc_ambi, opt.transition)
    bw = (opt.bw * 1.5 + 1.0).to_i32
    bw_long = (opt.bw_long * 1.5 + 1.0).to_i32
    bw_long = bw if bw_long < bw
    ep = MmExtra.new
    ep.dp_max2 = KSW_NEG_INF
    r.p = ep

    as1 = r.a_off
    cnt1 = r.cnt
    first_anchor = a[as1]
    last_anchor = a[as1 + cnt1 - 1]
    rs, qs = adjust_minier(mi, qseq0, first_anchor)
    re, qe = adjust_minier(mi, qseq0, last_anchor)

    seed_span = anchor_span(first_anchor)
    rs0 = u64_to_i32(first_anchor.x) + 1 - seed_span
    qs0 = u64_to_i32(first_anchor.y) + 1 - seed_span
    rs0 = 0 if rs0 < 0
    qs0 = 0 if qs0 < 0

    if qs > 0 && rs > 0
      lq = [qs, opt.max_gap].min
      qs0 = [qs0, qs - lq].min
      lr = lq
      lr += (lr * opt.a - opt.q) // opt.e if lr * opt.a > opt.q && opt.e > 0
      lr = [lr, opt.max_gap, rs].min
      rs0 = [rs0, rs - lr].min
    else
      rs0 = rs
      qs0 = qs
    end

    re0 = u64_to_i32(last_anchor.x) + 1
    qe0 = u64_to_i32(last_anchor.y) + 1
    if qe < qlen && re < ctg_len
      lq = [qlen - qe, opt.max_gap].min
      qe0 = [qe0, qe + lq].max
      lr = lq
      lr += (lr * opt.a - opt.q) // opt.e if lr * opt.a > opt.q && opt.e > 0
      lr = [lr, opt.max_gap, ctg_len - re].min
      re0 = [re0, re + lr].max
    else
      re0 = re
      qe0 = qe
    end
    qs0 = 0 if qs0 < 0
    rs0 = 0 if rs0 < 0
    qe0 = qlen if qe0 > qlen
    re0 = ctg_len if re0 > ctg_len

    rs1 = rs
    qs1 = qs
    if qs > qs0 && rs > rs0
      qpart = seq_window(qseq_all, qs0, qs)
      tpart = get_ref_seq(mi, rid, rs0, rs, false)
      seq_reverse(qpart)
      seq_reverse(tpart)
      ez = align_pair(opt, qpart.size, qpart, tpart.size, tpart, mat, bw,
        opt.end_bonus, r.split_inv? ? opt.zdrop_inv : opt.zdrop,
        KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT | KSW_EZ_REV_CIGAR)
      if ez.n_cigar > 0
        append_cigar(ep, ez.cigar)
        ep.dp_score += ez.max.to_i32
      end
      rs1 = rs - (ez.reach_end != 0 ? ez.mqe_t + 1 : ez.max_t + 1)
      qs1 = qs - (ez.reach_end != 0 ? qs - qs0 : ez.max_q + 1)
    end

    rs_cur = rs
    qs_cur = qs
    re1 = rs
    qe1 = qs
    i = 1
    while i < cnt1
      ai = a[as1 + i]
      if ((ai.y & (SEED_IGNORE | SEED_TANDEM)) != 0) && i != cnt1 - 1
        i += 1
        next
      end
      re_anchor, qe_anchor = adjust_minier(mi, qseq0, ai)
      re1 = re_anchor
      qe1 = qe_anchor
      if i == cnt1 - 1 || (ai.y & SEED_LONG_JOIN) != 0 ||
         (qe_anchor - qs_cur >= opt.min_ksw_len && re_anchor - rs_cur >= opt.min_ksw_len)
        bw1 = (ai.y & SEED_LONG_JOIN) != 0 ? [qe_anchor - qs_cur, re_anchor - rs_cur].max : bw_long
        qpart = seq_window(qseq_all, qs_cur, qe_anchor)
        tpart = get_ref_seq(mi, rid, rs_cur, re_anchor, false)
        if qpart.size > 0 && tpart.size > 0
          ez = align_pair(opt, qpart.size, qpart, tpart.size, tpart, mat, bw1,
            -1, opt.zdrop, KSW_EZ_APPROX_MAX)
          append_cigar(ep, ez.cigar) if ez.n_cigar > 0
          if ez.zdropped?
            ep.dp_score += ez.max.to_i32
            re1 = rs_cur + ez.max_t + 1
            qe1 = qs_cur + ez.max_q + 1
            break
          else
            ep.dp_score += ez.score
          end
        end
        rs_cur = re_anchor
        qs_cur = qe_anchor
      end
      i += 1
    end

    if qe1 < qe0 && re1 < re0
      qpart = seq_window(qseq_all, qe1, qe0)
      tpart = get_ref_seq(mi, rid, re1, re0, false)
      ez = align_pair(opt, qpart.size, qpart, tpart.size, tpart, mat, bw,
        opt.end_bonus, opt.zdrop, KSW_EZ_EXTZ_ONLY)
      if ez.n_cigar > 0
        append_cigar(ep, ez.cigar)
        ep.dp_score += ez.max.to_i32
      end
      re1 += ez.reach_end != 0 ? ez.mqe_t + 1 : ez.max_t + 1
      qe1 += ez.reach_end != 0 ? qe0 - qe1 : ez.max_q + 1
    end

    if ep.cigar.empty?
      r.p = nil
      return
    end

    r.rs = rs1
    r.re = re1
    if !rev
      r.qs = qs1
      r.qe = qe1
    else
      r.qs = qlen - qe1
      r.qe = qlen - qs1
    end

    full_qs = rev ? qlen - r.qe : r.qs
    full_qe = rev ? qlen - r.qs : r.qe
    return r.p = nil if full_qe <= full_qs || r.re <= r.rs
    q_aln = seq_window(qseq_all, full_qs, full_qe)
    t_aln = get_ref_seq(mi, rid, r.rs, r.re, false)
    update_extra(r, q_aln, t_aln, mat, opt.q, opt.e, (opt.flag & F_EQX) != 0)
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

    # Filter after alignment. filter_regs() compacts valid hits to the front
    # but does not resize the array; truncate so downstream code never sees
    # stale hits (see audit H-04).
    filter_regs(opt, qlen, n_regs_ref, regs)
    if n_regs_ref.value < regs.size
      regs.truncate(0, n_regs_ref.value)
    end

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
