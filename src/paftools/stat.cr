module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # stat — basic mapping statistics (PAF or SAM)
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_stat(args : Array(String)) : Int32
    gap_out_len : Int32? = nil; count_err = false
    rest = args.dup; positional = [] of String
    i = 0
    while i < rest.size
      case rest[i]
      when "-l"; i += 1; gap_out_len = rest[i].to_i
      when "-c"; count_err = true
      when "-h", "--help"
        STDERR.puts "Usage: paftools stat [-c] [-l INT] <in.sam>|<in.paf>"; return 0
      else positional << rest[i]
      end
      i += 1
    end
    if positional.empty?
      STDERR.puts "Usage: paftools stat [-c] [-l INT] <in.sam>|<in.paf>"; return 1
    end

    n_pri = 0_i64; n_2nd = 0_i64; n_seq = 0_i64; n_cigar_64k = 0_i64
    l_tot = 0_i64; l_cov = 0_i64; n_sub = 0_i64
    n_gap = [Array.new(6, 0_i64), Array.new(6, 0_i64)]
    last : String? = nil; last_qlen = 0
    regs = [] of {Int32, Int32}; lineno = 0

    open_in(positional[0]) do |io|
      io.each_line(chomp: true) do |line|
        lineno += 1; next if line.starts_with?('@')
        t = line.split('\t'); next if t.size < 2
        cigar : String? = nil; is_sam = false; is_rev = false; tname = ""
        rs = 0; mapq = 0; ori_qlen = 0; qs = 0; qe = 0
        nm : Int32? = nil; nn = 0

        if t.size > 4 && (t[4] == "+" || t[4] == "-" || t[4] == "*")
          next if t[4] == "*"
          if line.includes?("\ts2:i:")
            n_2nd += 1; next
          end
          nm = $1.to_i if line =~ /\tNM:i:(\d+)/
          nn = $1.to_i if line =~ /\tnn:i:(\d+)/
          if m = /\tcg:Z:(\S+)/.match(line)
            cigar = m[1]
          else
            STDERR.puts "WARNING: no CIGAR at line #{lineno}"; next
          end
          tname = t[5]; qs = t[2].to_i; qe = t[3].to_i
          is_rev = t[4] == "-"; rs = t[7].to_i
          mapq = t[11]?.try(&.to_i) || 0; ori_qlen = t[1].to_i
        else
          flag = t[1].to_i? || next
          next if (flag & 4) != 0 || t[2] == "*" || (t.size > 5 && t[5] == "*")
          if (flag & 0x100) != 0
            n_2nd += 1; next
          end
          nm = $1.to_i if line =~ /\tNM:i:(\d+)/
          nn = $1.to_i if line =~ /\tnn:i:(\d+)/
          cigar = t[5]; tname = t[2]; rs = t[3].to_i - 1; mapq = t[4].to_i
          is_sam = true; is_rev = (flag & 0x10) != 0
        end

        n_pri += 1
        if t[0] != last
          if last
            l_tot += last_qlen; l_cov += cov_len(regs)
          end
          regs = [] of {Int32, Int32}; n_seq += 1; last = t[0]
        end

        m_bases = 0; tl = 0; ql = 0; sclip = 0; clip = [0, 0]
        n_cop = 0; n_gapo = 0; n_gap_all = 0; l_match = 0

        cigar.try &.scan(/(\d+)([MIDNSHP=X])/) do |cm|
          l = cm[1].to_i; op = cm[2][0]; n_cop += 1
          case op
          when 'M', '=', 'X'; tl += l; ql += l; m_bases += l; l_match += l
          when 'I'
            ql += l
            bin = l < 50 ? 0 : l < 100 ? 1 : l < 300 ? 2 : l < 400 ? 3 : l < 1000 ? 4 : 5
            n_gap[0][bin] += 1
            puts [t[0], ql, is_rev ? '-' : '+', tname, rs + tl, 'I', l].join('\t') if gap_out_len.try { |g| l >= g }
            n_gapo += 1; n_gap_all += l
          when 'D'
            tl += l
            bin = l < 50 ? 0 : l < 100 ? 1 : l < 300 ? 2 : l < 400 ? 3 : l < 1000 ? 4 : 5
            n_gap[1][bin] += 1
            puts [t[0], ql, is_rev ? '-' : '+', tname, rs + tl, 'D', l].join('\t') if gap_out_len.try { |g| l >= g }
            n_gapo += 1; n_gap_all += l
          when 'N'; tl += l
          when 'S'; clip[m_bases == 0 ? 0 : 1] = l; sclip += l
          when 'H'; clip[m_bases == 0 ? 0 : 1] = l
          end
        end

        if nv = nm
          tmp = nv - n_gap_all - nn
          STDERR.puts "WARNING: NM < gaps at line #{lineno}" if tmp < 0 && nn == 0
          n_sub += [tmp, 0].max
        end
        n_cigar_64k += 1 if n_cop > 65535
        if is_sam
          qs = clip[is_rev ? 1 : 0]; qe = qs + ql; ori_qlen = clip[0] + ql + clip[1]
        end
        if count_err && (nv = nm)
          n_mm = [nv - n_gap_all, 0].max
          puts [t[0], ori_qlen, mapq, ori_qlen - (qe - qs), nv, l_match + n_gap_all, n_mm + n_gapo, l_match + n_gapo].join('\t')
        end
        regs << {qs, qe}; last_qlen = ori_qlen
      end
    end
    if regs.size > 0
      l_tot += last_qlen; l_cov += cov_len(regs)
    end

    unless gap_out_len || count_err
      puts "Number of mapped sequences: #{n_seq}"
      puts "Number of primary alignments: #{n_pri}"
      puts "Number of secondary alignments: #{n_2nd}"
      puts "Number of primary alignments with >65535 CIGAR operations: #{n_cigar_64k}"
      puts "Number of bases in mapped sequences: #{l_tot}"
      puts "Number of mapped bases: #{l_cov}"
      puts "Number of substitutions: #{n_sub}"
      [{"insertions", 0}, {"deletions", 1}].each do |(name, gi)|
        [[0, 50], [50, 100], [100, 300], [300, 400], [400, 1000], [-1, -1]].each_with_index do |(lo, hi), bi|
          range = hi < 0 ? "[1000,inf)" : "[#{lo},#{hi})"
          puts "Number of #{name} in #{range}: #{n_gap[gi][bi]}"
        end
      end
    end
    0
  end
end
