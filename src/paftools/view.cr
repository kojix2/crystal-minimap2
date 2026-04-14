module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # view — BLAST-like, MAF, or LASTZ-cigar output from PAF
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  private def self.pad(x, width, right = false) : String
    s = x.to_s
    return s if s.size >= width
    right ? s + " " * (width - s.size) : " " * (width - s.size) + s
  end

  # Append one cs operation to alignment string buffers.
  # elen[0/1] = cumulative ref/query lengths processed so far.
  private def self.update_aln(sref : String::Builder, sqry : String::Builder,
                              smid : String::Builder, type : Char, seq : String,
                              elen : Array(Int32)) : Nil
    case type
    when '=', ':'
      l = type == '=' ? seq.size : seq.to_i
      actual = type == '=' ? seq : seq # for ':' we don't have bases, use count
      if type == '='
        sref << seq; sqry << seq; smid << "|" * l
      else
        sref << " " * l; sqry << " " * l; smid << "|" * l
      end
      elen[0] += l; elen[1] += l
    when '*'
      sref << seq[0]; sqry << seq[1]; smid << " "
      elen[0] += 1; elen[1] += 1
    when '+'
      l = seq.size
      sref << "-" * l; sqry << seq; smid << " " * l
      elen[1] += l
    when '-'
      l = seq.size
      sref << seq; sqry << "-" * l; smid << " " * l
      elen[0] += l
    end
  end

  private def self.print_aln(rs : Int32, qs : Int32, strand : String,
                             slen : Array(Int32), elen : Array(Int32),
                             sref : String, sqry : String, smid : String) : Nil
    puts [pad("Ref+:", 5), pad(rs + slen[0] + 1, 10), sref, pad(rs + elen[0], 10, true)].join(" ")
    puts " " * 17 + smid
    if strand == "+"
      st = qs + slen[1] + 1; en = qs + elen[1]
    else
      st = qs - slen[1]; en = qs - elen[1] + 1
    end
    puts [pad("Qry#{strand}:", 5), pad(st, 10), sqry, pad(en, 10, true)].join(" ")
  end

  def self.cmd_view(args : Array(String)) : Int32
    line_len = 80; fmt = "aln"
    rest = [] of String
    i = 0
    while i < args.size
      case args[i]
      when "-f"; i += 1; fmt = args[i]
      when "-l"; i += 1; line_len = args[i].to_i
      when "-h", "--help"
        STDERR.puts "Usage: paftools view [-f aln|maf|lastz-cigar] [-l INT] <in.paf>"; return 0
      else rest << args[i]
      end
      i += 1
    end
    line_len = Int32::MAX if line_len == 0
    unless %w[aln maf lastz-cigar].includes?(fmt)
      STDERR.puts "Error: -f must be aln, maf, or lastz-cigar"; return 1
    end
    if rest.empty?
      STDERR.puts "Usage: paftools view [-f aln|maf|lastz-cigar] [-l INT] <in.paf>"; return 1
    end

    puts "##maf version=1\n" if fmt == "maf"
    lineno = 0

    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        lineno += 1
        t = line.split('\t', 12)
        next if t.size < 12

        case fmt
        when "lastz-cigar"
          cg = (m = /\tcg:Z:(\S+)/.match(line)) ? m[1] : nil
          unless cg
            STDERR.puts "WARNING: no cg tag at line #{lineno}"; next
          end
          score = (m = /\tAS:i:(\d+)/.match(line)) ? m[1] : "0"
          buf = ["cigar:", t[0], t[2], t[3], t[4], t[5], t[7], t[8], "+", score]
          cg.scan(/(\d+)([MIDNSHP=X])/) { |m| buf << m[2]; buf << m[1] }
          puts buf.join(" ")
        when "maf"
          cs = (m = /\tcs:Z:(\S+)/.match(line)) ? m[1] : nil
          unless cs
            STDERR.puts "WARNING: no cs tag at line #{lineno}"; next
          end
          sref = String::Builder.new; sqry = String::Builder.new
          smid = String::Builder.new; elen = [0, 0]
          cs.scan(/([:=*+\-])(\S+)/) do |m|
            update_aln(sref, sqry, smid, m[1][0], m[2], elen)
          end
          score = (m = /\tAS:i:(\d+)/.match(line)) ? m[1].to_i : 0
          len = [t[0].size, t[5].size].max
          ql = t[1].to_i
          qs, qe = t[4] == "+" ? {t[2].to_i, t[3].to_i} : {ql - t[3].to_i, ql - t[2].to_i}
          puts "a #{score}"
          puts [pad("s", 2), pad(t[5], len, true), pad(t[7], 10), pad(t[8].to_i - t[7].to_i, 10), "+", pad(t[6], 10, true), sref.to_s].join(" ")
          puts [pad("s", 2), pad(t[0], len, true), pad(qs, 10), pad(qe - qs, 10), t[4], pad(ql, 10, true), sqry.to_s].join(" ")
          puts ""
        else # aln (BLAST-like)
          cs = (m = /\tcs:Z:(\S+)/.match(line)) ? m[1] : nil
          unless cs
            STDERR.puts "WARNING: no cs tag at line #{lineno}"; next
          end
          # count stats
          n_mm = 0; n_oi = 0; n_od = 0; n_ei = 0; n_ed = 0
          cs.scan(/([:=*+\-])(\S+)/) do |m|
            case m[1][0]
            when '*'; n_mm += 1
            when '+'; n_oi += 1; n_ei += m[2].size
            when '-'; n_od += 1; n_ed += m[2].size
            end
          end
          bare = line.gsub(/\tc[sg]:Z:\S+/, "")
          puts ">#{bare}\tmm:i:#{n_mm}\toi:i:#{n_oi}\tei:i:#{n_ei}\tod:i:#{n_od}\ted:i:#{n_ed}"

          rs = t[7].to_i
          qs = t[4] == "+" ? t[2].to_i : t[3].to_i
          sref = ""; sqry = ""; smid = ""; sref_len = 0
          slen = [0, 0]; elen = [0, 0]; n_blocks = 0

          cs.scan(/([:=*+\-])(\S+)/) do |cm|
            op = cm[1][0]; seq = cm[2]
            rest_len = op == '*' ? 1 : (op == ':' ? seq.to_i : seq.size)
            start_pos = 0

            while rest_len > 0
              l_proc = [rest_len, line_len - sref_len].min
              # Build partial sequence for this block
              part = if op == '*'
                       seq
                     elsif op == ':'
                       " " * l_proc # short cs: no bases available, use spaces
                     else
                       seq[start_pos, l_proc]
                     end

              sb = String::Builder.new; qb = String::Builder.new; mb = String::Builder.new
              e2 = elen.dup
              case op
              when '=', ':'
                sb << part; qb << part; mb << "|" * l_proc
                e2[0] += l_proc; e2[1] += l_proc
              when '*'
                sb << seq[0]; qb << seq[1]; mb << " "
                e2[0] += 1; e2[1] += 1
              when '+'
                sb << "-" * l_proc; qb << part; mb << " " * l_proc; e2[1] += l_proc
              when '-'
                sb << part; qb << "-" * l_proc; mb << " " * l_proc; e2[0] += l_proc
              end
              sref += sb.to_s; sqry += qb.to_s; smid += mb.to_s
              sref_len += l_proc; elen = e2

              if sref_len >= line_len
                puts "" if n_blocks > 0
                print_aln(rs, qs, t[4], slen, elen, sref, sqry, smid)
                n_blocks += 1; sref = ""; sqry = ""; smid = ""; sref_len = 0
                slen = elen.dup
              end
              rest_len -= l_proc; start_pos += l_proc
            end
          end
          if sref_len > 0
            puts "" if n_blocks > 0
            print_aln(rs, qs, t[4], slen, elen, sref, sqry, smid)
          end
          puts "//"
        end
      end
    end
    0
  end
end
