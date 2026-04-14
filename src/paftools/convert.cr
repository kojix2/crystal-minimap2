module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # sam2paf — convert SAM to PAF
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_sam2paf(args : Array(String)) : Int32
    pri_only = false; pri_pri_only = false; allow_unmapped = false; long_cs = false
    rest = [] of String
    i = 0
    while i < args.size
      case args[i]
      when "-p"; pri_only = true
      when "-P"; pri_pri_only = true; pri_only = true
      when "-U"; allow_unmapped = true
      when "-L"; long_cs = true
      when "-h", "--help"
        STDERR.puts "Usage: paftools sam2paf [-pPUL] <in.sam>"; return 0
      else rest << args[i]
      end
      i += 1
    end
    if rest.empty?
      STDERR.puts "Usage: paftools sam2paf [-pPUL] <in.sam>"; return 1
    end

    ctg_len = Hash(String, Int32).new; lineno = 0

    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        lineno += 1; n_cigar = 0
        if line.starts_with?('@')
          if line.starts_with?("@SQ")
            name = (m = /\tSN:(\S+)/.match(line)) ? m[1] : nil
            l = (m = /\tLN:(\d+)/.match(line)) ? m[1].to_i : nil
            ctg_len[name] = l if name && l
          end
          next
        end
        t = line.split('\t', 11)
        flag = t[1].to_i
        if t.size < 10 || t[2] == "*" || (flag & 4) != 0 || t[5] == "*"
          if allow_unmapped
            qlen = t[9] == "*" ? 0 : t[9].size
            puts [t[0], qlen, 0, 0, "*", "*", 0, 0, 0, 0, 0, 0].join('\t')
          end
          next
        end
        next if pri_only && (flag & 0x100) != 0
        next if pri_pri_only && (flag & 0x900) != 0
        tlen = ctg_len[t[2]]?
        unless tlen
          STDERR.puts "Error at line #{lineno}: no length for contig #{t[2]}"; return 1
        end

        nm : Int32? = nil; md : String? = nil; cs_str : String? = nil
        if m = /\tNM:i:(\d+)/.match(line)
          nm = m[1].to_i
        end
        if m = /\tMD:Z:(\S+)/.match(line)
          md = m[1]
        end
        if m = /\tcs:Z:(\S+)/.match(line)
          cs_str = m[1]
        end
        md = nil if t[9] == "*"

        clip = [0, 0]; soft_clip = 0; i_ops = [0, 0]; d_ops = [0, 0]
        m_len = 0; n_span = 0; mm = 0; have_m = false; have_ext = false
        cigar_list = [] of {Int32, Char} # for MD reconstruction

        t[5].scan(/(\d+)([MIDNSHP=X])/) do |mat|
          l = mat[1].to_i; op = mat[2][0]; n_cigar += 1
          case op
          when 'M'; m_len += l; have_m = true
          when '='; m_len += l; have_ext = true
          when 'X'; m_len += l; mm += l; have_ext = true
          when 'I'; i_ops[0] += 1; i_ops[1] += l
          when 'D'; d_ops[0] += 1; d_ops[1] += l
          when 'N'; n_span += l
          when 'S'; clip[n_cigar == 1 ? 0 : 1] = l; soft_clip += l
          when 'H'; clip[n_cigar == 1 ? 0 : 1] = l
          end
          cigar_list << {l, op == '=' || op == 'X' ? 'M' : op} if md && op != 'H'
        end
        # Merge adjacent same-op entries in cigar_list
        merged = [] of {Int32, Char}
        cigar_list.each do |entry|
          if merged.empty? || merged.last[1] != entry[1]
            merged << entry
          else
            merged[-1] = {merged.last[0] + entry[0], entry[1]}
          end
        end

        ql_total = m_len + i_ops[1] + soft_clip
        tl_total = m_len + d_ops[1] + n_span
        ts = t[3].to_i - 1; te = ts + tl_total
        if n_cigar > 65535
          STDERR.puts "WARNING at line #{lineno}: #{n_cigar} CIGAR ops"
        end
        if te > tlen
          STDERR.puts "WARNING at line #{lineno}: alignment end > ref length; skipped"; next
        end
        if t[9] != "*" && t[9].size != ql_total
          STDERR.puts "WARNING at line #{lineno}: SEQ length inconsistent with CIGAR; skipped"; next
        end

        # Compute mm from CIGAR type
        if have_ext && !have_m
          nm_val = i_ops[1] + d_ops[1] + mm
          STDERR.puts "WARNING at line #{lineno}: NM differs" if nm && nm != nm_val
          nm = nm_val
        elsif nv = nm
          if nv < i_ops[1] + d_ops[1]
            STDERR.puts "WARNING at line #{lineno}: NM < total gaps"
            nm = i_ops[1] + d_ops[1]
          end
          mm = (nm || 0) - (i_ops[1] + d_ops[1])
        else
          STDERR.puts "WARNING at line #{lineno}: no NM, assuming 0 mismatches"
          mm = 0; nm = i_ops[1] + d_ops[1]
        end
        mlen_match = m_len - mm
        blen = m_len + i_ops[1] + d_ops[1]

        qlen = m_len + i_ops[1] + clip[0] + clip[1]
        qname = t[0]
        qname += "/1" if (flag & 1) != 0 && (flag & 0x40) != 0
        qname += "/2" if (flag & 1) != 0 && (flag & 0x80) != 0
        rev = (flag & 16) != 0
        qs = rev ? clip[1] : clip[0]; qe = qlen - (rev ? clip[0] : clip[1])

        # Build cs from MD if available and cs_str absent
        cs_parts = [] of String
        if md && cs_str.nil? && t[9] != "*"
          seq = t[9]
          k = 0; cx = 0; cy = 0; mx = 0; my = 0
          md.scan(/(\d+)|(\^[A-Za-z]+)|([A-Za-z])/) do |mm_m|
            if del_seq = mm_m[2]? # deletion ^XYZ
              len = del_seq.size - 1
              cs_parts << "-" << del_seq[1..]
              mx += len; cx += len; k += 1
            else
              ml = mm_m[1]? ? mm_m[1].to_i : 1
              while k < merged.size && merged[k][1] != 'D'
                cl = merged[k][0]; op = merged[k][1]
                if op == 'M'
                  if my + ml < cy + cl
                    if ml > 0
                      if mm_m[3]? # mismatch
                        cs_parts << "*" << mm_m[3] << seq[my..my]
                      elsif long_cs
                        cs_parts << "=" << seq[my, ml]
                      else
                        cs_parts << ":" << ml.to_s
                      end
                    end
                    mx += ml; my += ml; ml = 0; break
                  else
                    dl = cy + cl - my
                    cs_parts << (long_cs ? "=#{seq[my, dl]}" : ":#{dl}")
                    cx += cl; cy += cl; k += 1; mx += dl; my += dl; ml -= dl
                  end
                elsif op == 'I'
                  cs_parts << "+" << seq[cy, cl]
                  cy += cl; my += cl; k += 1
                elsif op == 'S'
                  cy += cl; my += cl; k += 1
                else
                  STDERR.puts "Warning: inconsistent MD at line #{lineno}"; break
                end
              end
            end
          end
        end

        type = (flag & 0x100) != 0 ? "S" : "P"
        tags = ["tp:A:#{type}"]
        if nv = nm
          tags << "NM:i:#{nv}"; tags << "mm:i:#{mm}"
        end
        tags << "gn:i:#{i_ops[1] + d_ops[1]}" << "go:i:#{i_ops[0] + d_ops[0]}"
        tags << "cg:Z:#{t[5].gsub(/\d+[SH]/, "")}"
        if cs_str
          tags << "cs:Z:#{cs_str}"
        elsif cs_parts.size > 0
          tags << "cs:Z:#{cs_parts.join}"
        end

        puts [qname, qlen, qs, qe, rev ? "-" : "+", t[2], tlen, ts, te,
              mlen_match, blen, t[4]].join('\t') + "\t" + tags.join('\t')
      end
    end
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # delta2paf — convert MUMmer delta to PAF
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_delta2paf(args : Array(String)) : Int32
    if args.empty? || args[0] == "-h" || args[0] == "--help"
      STDERR.puts "Usage: paftools delta2paf <in.delta>"; return args.empty? ? 1 : 0
    end
    rname = ""; qname = ""; rlen = 0; qlen = 0
    qs = 0; qe = 0; rs = 0; re_coord = 0; strand = 1; nm = 0
    cigar = [] of Int32 # packed: len<<4 | op (0=M,1=I,2=D)
    x = 0; y = 0; seen_gt = false

    open_in(args[0]) do |io|
      io.each_line(chomp: true) do |line|
        if m = /^>(\S+)\s+(\S+)\s+(\d+)\s+(\d+)/.match(line)
          rname = m[1]; qname = m[2]; rlen = m[3].to_i; qlen = m[4].to_i
          seen_gt = true; next
        end
        next unless seen_gt
        parts = line.split(' ')
        if parts.size == 7
          iparts = parts.map(&.to_i)
          r0 = iparts[0]; r1 = iparts[1]; q0 = iparts[2]; q1 = iparts[3]
          strand = ((r0 < r1 && q0 < q1) || (r0 > r1 && q0 > q1)) ? 1 : -1
          rs = (r0 < r1 ? r0 : r1) - 1; re_coord = r0 < r1 ? r1 : r0
          qs = (q0 < q1 ? q0 : q1) - 1; qe = q0 < q1 ? q1 : q0
          x = 0; y = 0; nm = iparts[4]; cigar = [] of Int32
        elsif parts.size == 1
          d = parts[0].to_i
          if d == 0
            raise "Inconsistent delta" if re_coord - rs - x != qe - qs - y
            cigar << ((re_coord - rs - x) << 4)
            blen = 0; cstr = [] of String
            cigar.each { |cig| blen += cig >> 4; cstr << "#{cig >> 4}#{"MID"[cig & 0xf]}" }
            puts [qname, qlen, qs, qe, strand > 0 ? "+" : "-",
                  rname, rlen, rs, re_coord, blen - nm, blen, 0,
                  "NM:i:#{nm}", "cg:Z:#{cstr.join}"].join('\t')
          elsif d > 0
            l = d - 1; x += l + 1; y += l
            cigar << (l << 4) if l > 0
            if !cigar.empty? && (cigar.last & 0xf) == 2
              cigar[-1] += 1 << 4
            else
              cigar << (1 << 4 | 2) # deletion
            end
          else
            l = -d - 1; x += l; y += l + 1
            cigar << (l << 4) if l > 0
            if !cigar.empty? && (cigar.last & 0xf) == 1
              cigar[-1] += 1 << 4
            else
              cigar << (1 << 4 | 1) # insertion
            end
          end
        end
      end
    end
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # longcs2seq — reconstruct ref or query sequence from long-cs PAF
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_longcs2seq(args : Array(String)) : Int32
    query_mode = false
    rest = [] of String
    i = 0
    while i < args.size
      case args[i]
      when "-q"          ; query_mode = true
      when "-h", "--help"; STDERR.puts "Usage: paftools longcs2seq [-q] <long-cs.paf>"; return 0
      else                 rest << args[i]
      end
      i += 1
    end
    if rest.empty?
      STDERR.puts "Usage: paftools longcs2seq [-q] <long-cs.paf>"; return 1
    end
    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('\t', 13)
        cs : String? = nil
        (12...t.size).each do |j|
          if m = /^cs:Z:(\S+)/.match(t[j])
            cs = m[1]; break
          end
        end
        next unless cs
        ts = String::Builder.new; qs_buf = String::Builder.new
        cs.scan(/([:=*+\-])([A-Za-z]+|\d+)/) do |mat|
          case mat[1][0]
          when '='; ts << mat[2]; qs_buf << mat[2]
          when '+'; qs_buf << mat[2].upcase
          when '-'; ts << mat[2].upcase
          when '*'; ts << mat[2][0].upcase; qs_buf << mat[2][1].upcase
          when ':'; raise "Long cs required (got short form ':'); rerun with --cs=long"
          end
        end
        if query_mode
          puts ">#{t[0]}_#{t[2]}_#{t[3]}"; puts qs_buf.to_s
        else
          puts ">#{t[5]}_#{t[7]}_#{t[8]}"; puts ts.to_s
        end
      end
    end
    0
  end
end
