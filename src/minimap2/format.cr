module Minimap2
  # ---------------------------------------------------------------------------
  # Output formatting — port of format.c
  # Supports PAF (primary), SAM (partial), cs/cigar strings.
  # ---------------------------------------------------------------------------

  BASES_STR = "ACGTN"
  BASES_LOWER = "acgtn"

  # Write CIGAR string (e.g. "10M3I5D") to a String::Builder
  private def self.write_cigar(sb : String::Builder, ep : MmExtra, is_eqx : Bool) : Nil
    ep.cigar.each do |c|
      op  = c & 0xf_u32
      len = c >> 4
      sb << len
      sb << CIGAR_STR[op.to_i]?
    end
  end

  # Compute CIGAR length in query / target bases
  private def self.cigar_qlen_tlen(ep : MmExtra) : {Int32, Int32}
    q = 0; t = 0
    ep.cigar.each do |c|
      op  = (c & 0xf).to_i32
      len = (c >> 4).to_i32
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

  # Write cs string (short form with :N for matches).
  private def self.write_cs(sb : String::Builder, tseq : Array(UInt8), qseq : Array(UInt8),
                              ep : MmExtra, long_form : Bool) : Nil
    q_off = 0; t_off = 0
    ep.cigar.each do |c|
      op  = (c & 0xf).to_i32
      len = (c >> 4).to_i32
      case op
      when CIGAR_MATCH, CIGAR_EQ_MATCH, CIGAR_X_MISMATCH
        l_tmp = 0
        len.times do |j|
          if qseq[q_off + j] != tseq[t_off + j]
            if l_tmp > 0
              if long_form
                sb << "="; sb << BASES_STR[0, 1] # placeholder
              else
                sb << ":"; sb << l_tmp
              end
              l_tmp = 0
            end
            sb << "*"
            sb << BASES_LOWER[tseq[t_off + j]]?
            sb << BASES_LOWER[qseq[q_off + j]]?
          else
            l_tmp += 1
          end
        end
        if l_tmp > 0
          if long_form
            sb << "="; l_tmp.times { |j| sb << BASES_STR[qseq[q_off + j + (len - l_tmp)]]? }
          else
            sb << ":"; sb << l_tmp
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
        d = tseq[t_off]? || 4_u8; a = tseq[t_off + len - 1]? || 4_u8
        sb << "~"; sb << BASES_LOWER[d]? ; sb << BASES_LOWER[tseq[t_off + 1]? || 4_u8]?
        sb << len; sb << BASES_LOWER[tseq[t_off + len - 2]? || 4_u8]?; sb << BASES_LOWER[a]?
        t_off += len
      end
    end
  end

  # Write one PAF record.
  # Mirrors mm_write_paf4().
  def self.write_paf(io : IO, mi : MmIdx, t : BSeq1, r : MmReg1,
                      opt_flag : Int64, rep_len : Int32 = 0) : Nil
    ep = r.p

    # Column 1-12: standard PAF fields
    ctg_name = r.rid >= 0 && r.rid < mi.seq.size ? mi.seq[r.rid].name : "*"
    ctg_len  = r.rid >= 0 && r.rid < mi.seq.size ? mi.seq[r.rid].len.to_i32 : 0

    mapq = r.mapq.clamp(0_u32, 60_u32)

    io.print t.name
    io.print "\t"; io.print t.l_seq
    io.print "\t"; io.print r.qs
    io.print "\t"; io.print r.qe
    io.print "\t"; io.print r.rev ? '-' : '+'
    io.print "\t"; io.print ctg_name
    io.print "\t"; io.print ctg_len
    io.print "\t"; io.print r.rs
    io.print "\t"; io.print r.re
    io.print "\t"; io.print r.mlen
    io.print "\t"; io.print r.blen
    io.print "\t"; io.print mapq

    # Optional fields
    if r.parent >= 0 && r.parent != r.id
      io.print "\ttp:A:S"
    else
      io.print "\ttp:A:P"
    end
    io.print "\tcm:i:"; io.print r.cnt
    io.print "\ts1:i:"; io.print r.score
    if ep
      io.print "\ts2:i:"; io.print ep.dp_max2
    end
    if r.div >= 0.0_f32
      io.printf("\tdv:f:%.4f", r.div)
    end
    if ep && (opt_flag & F_OUT_MD) == 0 && (opt_flag & F_CIGAR) != 0
      # Write CIGAR
      io.print "\tcg:Z:"
      sb = String::Builder.new
      write_cigar(sb, ep, (opt_flag & F_EQX) != 0)
      io.print sb.to_s
    end
    if (opt_flag & F_OUT_MD) != 0 && ep
      # MD tag skipped for now
    end
    io.print "\n"
  end

  # Write SAM header.
  def self.write_sam_hdr(io : IO, mi : MmIdx?, version : String = LIB_VERSION) : Nil
    io.puts "@HD\tVN:1.6\tSO:unsorted\tGO:query"
    if mi
      mi.seq.each do |s|
        io.puts "@SQ\tSN:#{s.name}\tLN:#{s.len}"
      end
    end
    io.puts "@PG\tID:minimap2\tPN:minimap2\tVN:#{version}"
  end

  # Write one SAM record.
  def self.write_sam(io : IO, mi : MmIdx, t : BSeq1, r : MmReg1, n_regs : Int32,
                     all_regs : Array(MmReg1), opt_flag : Int64) : Nil
    ep = r.p
    ctg_name = r.rid >= 0 && r.rid < mi.seq.size ? mi.seq[r.rid].name : "*"
    mapq = r.mapq.clamp(0_u32, 60_u32)
    flag = 0_u32

    flag |= 0x10_u32 if r.rev            # reverse complement
    flag |= 0x100_u32 unless r.sam_pri   # not primary
    flag |= 0x800_u32 if r.is_alt        # supplementary

    io.print t.name; io.print "\t"; io.print flag
    io.print "\t"; io.print ctg_name
    io.print "\t"; io.print r.rs + 1  # 1-based
    io.print "\t"; io.print mapq
    io.print "\t"

    if ep
      sb = String::Builder.new
      # soft-clip at start
      if r.qs > 0
        sb << r.qs; sb << "S"
      end
      write_cigar(sb, ep, false)
      if t.l_seq - r.qe > 0
        sb << (t.l_seq - r.qe); sb << "S"
      end
      io.print sb.to_s
    else
      io.print "*"
    end

    io.print "\t*\t0\t0\t"  # RNEXT, PNEXT, TLEN

    # SEQ
    if (opt_flag & F_NO_QUAL) == 0
      io.print t.seq
    else
      io.print "*"
    end

    io.print "\t"
    # QUAL
    if t.qual
      if r.rev
        io.print t.qual.not_nil!.reverse
      else
        io.print t.qual
      end
    else
      io.print "*"
    end

    io.print "\n"
  end
end
