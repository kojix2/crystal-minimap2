module Minimap2
  # ---------------------------------------------------------------------------
  # Output formatting — port of format.c
  # Supports PAF (primary), SAM, cs/cigar strings, MD tag.
  # ---------------------------------------------------------------------------

  BASES_STR   = "ACGTN"
  BASES_LOWER = "acgtn"

  # Write CIGAR string (e.g. "10M3I5D") to a String::Builder
  private def self.write_cigar(sb : String::Builder, ep : MmExtra, is_eqx : Bool) : Nil
    ep.cigar.each do |entry|
      op = entry & 0xf_u32
      len = entry >> 4
      sb << len
      sb << CIGAR_STR[op.to_i]?
    end
  end

  # Compute CIGAR length in query / target bases
  private def self.cigar_qlen_tlen(ep : MmExtra) : {Int32, Int32}
    q = 0; t = 0
    ep.cigar.each do |entry|
      op = (entry & 0xf).to_i32
      len = (entry >> 4).to_i32
      case op
      when CIGAR_MATCH, CIGAR_EQ_MATCH, CIGAR_X_MISMATCH
        q += len; t += len
      when CIGAR_INS, CIGAR_SOFTCLIP
        q += len
      when CIGAR_DEL, CIGAR_N_SKIP
        t += len
      end
    end
    {q, t}
  end

  # Write cs string (short or long form).
  # Short: :N *tc +acgt -acgt ~gt..ac
  # Long:  =ACGT *tc +acgt -acgt ~gt..ac
  private def self.write_cs(sb : String::Builder, tseq : Array(UInt8), qseq : Array(UInt8),
                            ep : MmExtra, long_form : Bool) : Nil
    q_off = 0; t_off = 0
    ep.cigar.each do |entry|
      op = (entry & 0xf).to_i32
      len = (entry >> 4).to_i32
      case op
      when CIGAR_MATCH, CIGAR_EQ_MATCH, CIGAR_X_MISMATCH
        # Track matching and mismatching runs
        match_start = 0
        j = 0
        while j < len
          if qseq[q_off + j] == tseq[t_off + j]
            # Start or extend a match run
            match_start = j if j == 0 || qseq[q_off + j - 1] != tseq[t_off + j - 1]
          else
            # End any preceding match run
            match_len = j - match_start
            if qseq[q_off + j - 1] == tseq[t_off + j - 1]
              # Previous was match
              if long_form
                sb << "="
                match_len.times { |k| sb << BASES_STR[qseq[q_off + match_start + k]]? }
              else
                sb << ":"; sb << match_len
              end
            end
            # Write substitution
            sb << "*"
            sb << BASES_LOWER[tseq[t_off + j]]?
            sb << BASES_LOWER[qseq[q_off + j]]?
            match_start = j + 1
          end
          j += 1
        end
        # Handle trailing match
        final_match_len = len - match_start
        if final_match_len > 0
          if long_form
            sb << "="
            final_match_len.times { |k| sb << BASES_STR[qseq[q_off + match_start + k]]? }
          else
            sb << ":"; sb << final_match_len
          end
        end
        q_off += len; t_off += len
      when CIGAR_INS
        sb << "+"
        len.times { |j| sb << BASES_LOWER[qseq[q_off + j]]? }
        q_off += len
      when CIGAR_DEL
        sb << "-"
        len.times { |j| sb << BASES_LOWER[tseq[t_off + j]]? }
        t_off += len
      when CIGAR_N_SKIP
        d = tseq[t_off]? || 4_u8
        a = tseq[t_off + len - 1]? || 4_u8
        sb << "~"
        sb << BASES_LOWER[d]?
        sb << BASES_LOWER[tseq[t_off + 1]? || 4_u8]?
        sb << len
        sb << BASES_LOWER[tseq[t_off + len - 2]? || 4_u8]?
        sb << BASES_LOWER[a]?
        t_off += len
      end
    end
  end

  # Write MD tag from CIGAR and sequences.
  # Mirrors write_md_core() from format.c.
  private def self.write_md(sb : String::Builder, tseq : Array(UInt8), qseq : Array(UInt8),
                            ep : MmExtra) : Nil
    q_off = 0; t_off = 0; l_md = 0
    ep.cigar.each do |entry|
      op = (entry & 0xf).to_i32
      len = (entry >> 4).to_i32
      case op
      when CIGAR_MATCH, CIGAR_EQ_MATCH, CIGAR_X_MISMATCH
        len.times do |j|
          if qseq[q_off + j] == tseq[t_off + j]
            l_md += 1
          else
            sb << l_md
            sb << BASES_STR[tseq[t_off + j]]?
            l_md = 0
          end
        end
        q_off += len; t_off += len
      when CIGAR_INS
        q_off += len # insertions don't affect MD
      when CIGAR_DEL
        sb << l_md; l_md = 0
        sb << "^"
        len.times { |j| sb << BASES_STR[tseq[t_off + j]]? }
        t_off += len
      when CIGAR_N_SKIP
        t_off += len # skips don't affect MD
      end
    end
    sb << l_md if l_md > 0
  end

  # Compute edit distance (NM tag) from CIGAR.
  private def self.compute_nm(ep : MmExtra) : Int32
    n_mis = 0
    q_off = 0; t_off = 0
    ep.cigar.each do |entry|
      op = (entry & 0xf).to_i32
      len = (entry >> 4).to_i32
      case op
      when CIGAR_MATCH, CIGAR_EQ_MATCH, CIGAR_X_MISMATCH
        len.times { |j| n_mis += 1 } # will be corrected by counting matches
      when CIGAR_INS, CIGAR_DEL
        n_mis += len
      end
    end
    # NM = blen - mlen + n_ambi
    n_mis = 0
    ep.cigar.each do |entry|
      op = (entry & 0xf).to_i32
      len = (entry >> 4).to_i32
      if op == CIGAR_INS || op == CIGAR_DEL
        n_mis += len
      end
    end
    n_mis
  end

  # Write one PAF record.
  # Mirrors mm_write_paf4().
  def self.write_paf(io : IO, mi : MmIdx, t : BSeq1, r : MmReg1,
                     opt_flag : Int64, rep_len : Int32 = 0,
                     tseq : Array(UInt8)? = nil, qseq : Array(UInt8)? = nil) : Nil
    ep = r.p

    # Column 1-12: standard PAF fields
    ctg_name = r.rid >= 0 && r.rid < mi.seq.size ? mi.seq[r.rid].name : "*"
    ctg_len = r.rid >= 0 && r.rid < mi.seq.size ? mi.seq[r.rid].len.to_i32 : 0

    mapq = r.mapq.clamp(0_u32, 60_u32)

    io.print t.name
    io.print "\t"; io.print t.l_seq
    io.print "\t"; io.print r.qs
    io.print "\t"; io.print r.qe
    io.print "\t"; io.print r.rev? ? '-' : '+'
    io.print "\t"; io.print ctg_name
    io.print "\t"; io.print ctg_len
    if (opt_flag & F_QSTRAND) != 0 && r.rev?
      io.print "\t"; io.print ctg_len - r.re
      io.print "\t"; io.print ctg_len - r.rs
    else
      io.print "\t"; io.print r.rs
      io.print "\t"; io.print r.re
    end
    io.print "\t"; io.print r.mlen
    io.print "\t"; io.print r.blen
    io.print "\t"; io.print mapq

    # Optional fields — order matches C minimap2 write_tags()
    # NM, ms, AS, nn tags (when CIGAR/alignment exists)
    if ep
      nm = r.blen - r.mlen + ep.n_ambi.to_i32
      io.print "\tNM:i:"; io.print nm
      io.print "\tms:i:"; io.print ep.dp_max0
      io.print "\tAS:i:"; io.print ep.dp_score
      io.print "\tnn:i:"; io.print ep.n_ambi
    end

    # ts tag (transcript strand)
    if ep && ep.trans_strand > 0
      io.print "\tts:A:"
      case ep.trans_strand
      when 1 then io.print "+"
      when 2 then io.print "-"
      else        io.print "?"
      end
    end

    # tp tag — primary/secondary inversion distinction (mirrors write_tags())
    if r.inv?
      io.print "\ttp:A:"
      io.print (r.parent >= 0 && r.parent != r.id) ? "i" : "I"
    elsif r.parent >= 0 && r.parent != r.id
      io.print "\ttp:A:S"
    else
      io.print "\ttp:A:P"
    end

    # cm tag
    io.print "\tcm:i:"; io.print r.cnt
    # s1 tag
    io.print "\ts1:i:"; io.print r.score
    # s2 tag (primary alignments only)
    if r.parent == r.id
      io.print "\ts2:i:"; io.print r.subsc
    end

    # de/dv tag (divergence)
    if ep
      # de:f: — event-level divergence from CIGAR
      n_gap = 0
      n_gapo = 0
      ep.cigar.each do |entry|
        op = (entry & 0xf).to_i32
        len = (entry >> 4).to_i32
        if op == CIGAR_INS || op == CIGAR_DEL
          n_gap += len
          n_gapo += 1
        end
      end
      denom = r.blen + ep.n_ambi.to_i32 - n_gap + n_gapo
      if denom > 0
        div = 1.0 - r.mlen.to_f / denom
        io.printf("\tde:f:%.4f", div)
      else
        io.print "\tde:f:0"
      end
    elsif r.div >= 0.0_f32 && r.div <= 1.0_f32
      # dv:f: — approximate divergence from minimizer matching
      if r.div == 0.0_f32
        io.print "\tdv:f:0"
      else
        io.printf("\tdv:f:%.4f", r.div)
      end
    end

    # zd tag (split alignment)
    if r.split > 0
      io.print "\tzd:i:"; io.print r.split
    end

    # rl tag (repetitive length)
    if rep_len >= 0
      io.print "\trl:i:"; io.print rep_len
    end

    # cg tag (CIGAR) — gated only on F_OUT_CG, like C minimap2
    if ep && (opt_flag & F_OUT_CG) != 0
      io.print "\tcg:Z:"
      sb = String::Builder.new
      write_cigar(sb, ep, (opt_flag & F_EQX) != 0)
      io.print sb.to_s
    end

    # cs tag
    if ep && tseq && qseq && ((opt_flag & F_OUT_CS) != 0 || (opt_flag & F_OUT_CS_LONG) != 0)
      io.print "\tcs:Z:"
      sb = String::Builder.new
      write_cs(sb, tseq, qseq, ep, (opt_flag & F_OUT_CS_LONG) != 0)
      io.print sb.to_s
    end

    # MD tag
    if ep && tseq && qseq && (opt_flag & F_OUT_MD) != 0
      io.print "\tMD:Z:"
      sb = String::Builder.new
      write_md(sb, tseq, qseq, ep)
      io.print sb.to_s
    end

    if (opt_flag & F_COPY_COMMENT) != 0 && (comment = t.comment)
      io.print "\t"; io.print comment
    end

    io.print "\n"
  end

  # Write SAM header. If `rg_line` (from -R) is given it is emitted verbatim
  # and the trailing `ID:...` value is reused as the RG tag for each record.
  def self.write_sam_hdr(io : IO, mi : MmIdx?, version : String = LIB_VERSION,
                         rg_line : String? = nil) : Nil
    io.puts "@HD\tVN:1.6\tSO:unsorted\tGO:query"
    if mi
      mi.seq.each do |seq_rec|
        io.puts "@SQ\tSN:#{seq_rec.name}\tLN:#{seq_rec.len}"
      end
    end
    if rg_line
      normalized_rg = rg_line.gsub("\\t", "\t")
      io.puts normalized_rg.starts_with?("@RG") ? normalized_rg : "@RG\t#{normalized_rg}"
    end
    io.puts "@PG\tID:minimap2\tPN:minimap2\tVN:#{version}"
  end

  # Extract the RG ID from a -R header line ("@RG\tID:foo\tSM:bar" → "foo").
  def self.rg_id(rg_line : String?) : String?
    return nil unless rg_line
    rg_line.gsub("\\t", "\t").split('\t').each do |field|
      if field.starts_with?("ID:")
        return field[3..]
      end
    end
    nil
  end

  # Write one SAM record.
  def self.write_sam(io : IO, mi : MmIdx, t : BSeq1, r : MmReg1, n_regs : Int32,
                     all_regs : Array(MmReg1), opt_flag : Int64,
                     tseq : Array(UInt8)? = nil, qseq : Array(UInt8)? = nil,
                     rg_id : String? = nil) : Nil
    ep = r.p
    ctg_name = r.rid >= 0 && r.rid < mi.seq.size ? mi.seq[r.rid].name : "*"
    mapq = r.mapq.clamp(0_u32, 60_u32)
    flag = 0_u32

    flag |= 0x10_u32 if r.rev?          # reverse complement
    # secondary (0x100): parent != id
    # supplementary (0x800): parent == id but not sam_pri
    # (mirrors mm_write_sam3() in format.c; is_alt is unrelated to SAM flags)
    if r.parent >= 0 && r.parent != r.id
      flag |= 0x100_u32
    elsif !r.sam_pri?
      flag |= 0x800_u32
    end

    io.print t.name; io.print "\t"; io.print flag
    io.print "\t"; io.print ctg_name
    io.print "\t"; io.print r.rs + 1 # 1-based
    io.print "\t"; io.print mapq
    io.print "\t"

    if ep
      sb = String::Builder.new
      # soft-clip: for reverse-strand alignments the 5'/3' clips swap
      # (mirrors write_sam_cigar() in format.c).
      clip5 = r.rev? ? t.l_seq - r.qe : r.qs
      clip3 = r.rev? ? r.qs : t.l_seq - r.qe
      if clip5 > 0
        sb << clip5; sb << "S"
      end
      write_cigar(sb, ep, false)
      if clip3 > 0
        sb << clip3; sb << "S"
      end
      io.print sb.to_s
    else
      io.print "*"
    end

    io.print "\t*\t0\t0\t" # RNEXT, PNEXT, TLEN (single-segment only)

    # SEQ: for reverse-strand alignments output the reverse-complement.
    if r.rev?
      io.print t.seq.reverse.tr("ACGTUacgtuNn", "TGCAAtgcaaNn")
    else
      io.print t.seq
    end

    io.print "\t"
    # QUAL: for reverse-strand alignments output the reversed quality string.
    if t.qual && (opt_flag & F_NO_QUAL) == 0
      if r.rev?
        io.print t.qual.try(&.reverse)
      else
        io.print t.qual
      end
    else
      io.print "*"
    end

    # RG tag (read group) — placed before other tags like C minimap2
    if rg_id
      io.print "\tRG:Z:"; io.print rg_id
    end

    # Optional tags
    if ep
      # NM tag
      nm = r.blen - r.mlen + ep.n_ambi.to_i32
      io.print "\tNM:i:"; io.print nm

      # ms tag
      io.print "\tms:i:"; io.print ep.dp_max0

      # AS tag
      io.print "\tAS:i:"; io.print ep.dp_score

      # nn tag
      io.print "\tnn:i:"; io.print ep.n_ambi

      # ts tag (transcript strand)
      if ep.trans_strand > 0
        io.print "\tts:A:"
        case ep.trans_strand
        when 1 then io.print "+"
        when 2 then io.print "-"
        else        io.print "?"
        end
      end

      # cs tag
      if tseq && qseq && ((opt_flag & F_OUT_CS) != 0 || (opt_flag & F_OUT_CS_LONG) != 0)
        io.print "\tcs:Z:"
        sb = String::Builder.new
        write_cs(sb, tseq, qseq, ep, (opt_flag & F_OUT_CS_LONG) != 0)
        io.print sb.to_s
      end

      # MD tag
      if tseq && qseq && (opt_flag & F_OUT_MD) != 0
        io.print "\tMD:Z:"
        sb = String::Builder.new
        write_md(sb, tseq, qseq, ep)
        io.print sb.to_s
      end
    end

    if (opt_flag & F_COPY_COMMENT) != 0 && (comment = t.comment)
      io.print "\t"; io.print comment
    end

    io.print "\n"
  end
end
