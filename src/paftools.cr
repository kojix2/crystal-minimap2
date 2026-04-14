require "option_parser"
require "compress/gzip"

# paftools.cr — Crystal port of paftools.js (minimap2 companion utilities)
module Paftools
  VERSION = "2.30 (Crystal port)"

  # ── I/O ──────────────────────────────────────────────────────────────────

  private def self.open_in(fn : String, &)
    if fn == "-"
      yield STDIN
    elsif fn.ends_with?(".gz") || fn.ends_with?(".bgz")
      File.open(fn, "rb") { |f| Compress::Gzip::Reader.open(f) { |gz| yield gz } }
    else
      File.open(fn) { |f| yield f }
    end
  end

  # ── Interval operations ───────────────────────────────────────────────────
  # Sorted {start,end} pairs with a link-back index for O(n) overlap queries
  # (mirrors Interval.sort / merge / index_end / find_ovlp in paftools.js).

  def self.intv_sort(a : Array({Int32, Int32}))
    a.sort_by! { |iv| {iv[0], iv[1]} }
  end

  def self.intv_merge(a : Array({Int32, Int32}))
    return if a.size < 2
    k = 0
    (1...a.size).each do |i|
      if a[k][1] >= a[i][0]
        a[k] = {a[k][0], [a[k][1], a[i][1]].max}
      else
        k += 1; a[k] = a[i]
      end
    end
    a.delete_at(k + 1, a.size - k - 1) if a.size > k + 1
  end

  def self.intv_build(a : Array({Int32, Int32})) : Array(Int32)
    n = a.size; return [] of Int32 if n == 0
    idx = Array(Int32).new(n, 0)
    k = 0; k_en = a[0][1]
    (1...n).each do |i|
      if k_en <= a[i][0]
        k += 1
        while k < i
          break if a[k][1] > a[i][0]; k += 1
        end
        k_en = a[k][1]
      end
      idx[i] = k
    end
    idx
  end

  # Returns indices into `a` of all intervals overlapping [st, en)
  def self.intv_ovlp_idx(a : Array({Int32, Int32}), idx : Array(Int32),
                         st : Int32, en : Int32) : Array(Int32)
    return [] of Int32 if a.empty? || st >= en
    left = -1; right = a.size
    while right - left > 1
      mid = left + ((right - left) >> 1)
      if a[mid][0] > st
        right = mid
      elsif a[mid][0] < st
        left = mid
      else
        left = mid; break
      end
    end
    k = left < 0 ? 0 : idx[left]
    res = [] of Int32
    i = k
    while i < a.size
      break if a[i][0] >= en
      res << i if st < a[i][1]
      i += 1
    end
    res
  end

  def self.intv_ovlp(a : Array({Int32, Int32}), idx : Array(Int32),
                     st : Int32, en : Int32) : Array({Int32, Int32})
    intv_ovlp_idx(a, idx, st, en).map { |i| a[i] }
  end

  # ── Merge overlapping regions and return covered length ───────────────────

  private def self.cov_len(regs : Array({Int32, Int32})) : Int64
    return 0_i64 if regs.empty?
    s = regs.sort_by { |r| r[0] }
    st = s[0][0]; en = s[0][1]; l = 0_i64
    (1...s.size).each do |i|
      if s[i][0] < en
        en = s[i][1] > en ? s[i][1] : en
      else
        l += en - st; st = s[i][0]; en = s[i][1]
      end
    end
    l + en - st
  end

  # ── Tiny FASTA reader (for cmd_call -f) ──────────────────────────────────

  private def self.read_fasta(fn : String) : {Hash(String, String), Array({String, Int32})}
    h = Hash(String, String).new; lens = [] of {String, Int32}
    open_in(fn) do |io|
      name = ""; seq = String::Builder.new
      io.each_line(chomp: true) do |line|
        if line.starts_with?('>')
          unless name.empty?
            s = seq.to_s; h[name] = s; lens << {name, s.size}
          end
          name = line[1..].split(' ', 2)[0]; seq = String::Builder.new
        else
          seq << line
        end
      end
      unless name.empty?
        s = seq.to_s; h[name] = s; lens << {name, s.size}
      end
    end
    {h, lens}
  end

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
        rs = 0; mapq = 0; aqlen = 0; ori_qlen = 0; qs = 0; qe = 0
        atlen : Int32? = nil; nm : Int32? = nil; nn = 0

        if t.size > 4 && (t[4] == "+" || t[4] == "-" || t[4] == "*")
          next if t[4] == "*"
          if line.includes?("\ts2:i:")
            n_2nd += 1; next
          end
          nm = $1.to_i if (line =~ /\tNM:i:(\d+)/)
          nn = $1.to_i if (line =~ /\tnn:i:(\d+)/)
          if (m = /\tcg:Z:(\S+)/.match(line))
            cigar = m[1]
          else
            STDERR.puts "WARNING: no CIGAR at line #{lineno}"; next
          end
          tname = t[5]; qs = t[2].to_i; qe = t[3].to_i; aqlen = qe - qs
          is_rev = t[4] == "-"; rs = t[7].to_i; atlen = t[8].to_i - rs
          mapq = t[11]?.try(&.to_i) || 0; ori_qlen = t[1].to_i
        else
          flag = t[1].to_i? || next
          next if (flag & 4) != 0 || t[2] == "*" || (t.size > 5 && t[5] == "*")
          if (flag & 0x100) != 0
            n_2nd += 1; next
          end
          nm = $1.to_i if (line =~ /\tNM:i:(\d+)/)
          nn = $1.to_i if (line =~ /\tnn:i:(\d+)/)
          cigar = t[5]; tname = t[2]; rs = t[3].to_i - 1; mapq = t[4].to_i
          aqlen = t[9].size; is_sam = true; is_rev = (flag & 0x10) != 0
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

        if (nv = nm)
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

        nn = 0; nm : Int32? = nil; md : String? = nil; cs_str : String? = nil
        if (m = /\tNM:i:(\d+)/.match(line))
          nm = m[1].to_i
        end
        if (m = /\tnn:i:(\d+)/.match(line))
          nn = m[1].to_i
        end
        if (m = /\tMD:Z:(\S+)/.match(line))
          md = m[1]
        end
        if (m = /\tcs:Z:(\S+)/.match(line))
          cs_str = m[1]
        end
        md = nil if t[9] == "*"

        clip = [0, 0]; soft_clip = 0; i_ops = [0, 0]; d_ops = [0, 0]
        m_len = 0; n_span = 0; mm = 0; have_m = false; have_ext = false
        cigar_list = [] of {Int32, Char} # for MD reconstruction

        t[5].scan(/(\d+)([MIDNSHP=X])/) do |cm|
          l = cm[1].to_i; op = cm[2][0]; n_cigar += 1
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
        elsif (nv = nm)
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
            if (del_seq = mm_m[2]?) # deletion ^XYZ
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
        if (nv = nm)
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
        if (m = /^>(\S+)\s+(\S+)\s+(\d+)\s+(\d+)/.match(line))
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
            cigar.each { |c| blen += c >> 4; cstr << "#{c >> 4}#{"MID"[c & 0xf]}" }
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
          if (m = /^cs:Z:(\S+)/.match(t[j]))
            cs = m[1]; break
          end
        end
        next unless cs
        ts = String::Builder.new; qs_buf = String::Builder.new
        cs.scan(/([:=*+\-])([A-Za-z]+|\d+)/) do |m|
          case m[1][0]
          when '='; ts << m[2]; qs_buf << m[2]
          when '+'; qs_buf << m[2].upcase
          when '-'; ts << m[2].upcase
          when '*'; ts << m[2][0].upcase; qs_buf << m[2][1].upcase
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

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # liftover — lift BED coordinates via PAF CIGAR
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_liftover(args : Array(String)) : Int32
    to_merge = false; min_mapq = 5; min_len = 50000; max_div = 2.0
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-m"; to_merge = true
      when "-q"; i += 1; min_mapq = args[i].to_i
      when "-l"; i += 1; min_len = args[i].to_i
      when "-d"; i += 1; max_div = args[i].to_f
      when "-h", "--help"
        STDERR.puts "Usage: paftools liftover [-m] [-q INT] [-l INT] [-d FLOAT] <aln.paf> <query.bed>"
        return 0
      else rest << args[i]
      end
      i += 1
    end
    if rest.size < 2
      STDERR.puts "Usage: paftools liftover [-m] [-q INT] [-l INT] [-d FLOAT] <aln.paf> <query.bed>"
      return 1
    end

    # Read BED into per-chrom sorted (and optionally merged) interval lists
    bed = Hash(String, Array({Int32, Int32})).new
    bed_idx = Hash(String, Array(Int32)).new
    open_in(rest[1]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('\t'); next if t.size < 3
        chr = t[0]
        bed[chr] ||= [] of {Int32, Int32}
        bed[chr] << {t[1].to_i, t[2].to_i}
      end
    end
    bed.each do |chr, ivs|
      intv_sort(ivs)
      intv_merge(ivs) if to_merge
      bed_idx[chr] = intv_build(ivs)
    end

    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('\t')
        next if t.size < 12
        ivs = bed[t[0]]?; next unless ivs
        tp : String? = nil; cg : String? = nil
        t[12..].each do |tag|
          if (m = /^(\S\S):[AZif]:(\S+)/.match(tag))
            tp = m[2] if m[1] == "tp"
            cg = m[2] if m[1] == "cg"
          end
        end
        next unless tp == "P" || tp == "I"
        next unless cg
        (1..3).each { |j| t[j] = t[j].to_i.to_s }
        (6..11).each { |j| t[j] = t[j].to_i.to_s }
        next if t[11].to_i < min_mapq || t[10].to_i < min_len

        regs = intv_ovlp(ivs, bed_idx[t[0]], t[2].to_i, t[3].to_i)
        next if regs.empty?

        if max_div >= 0.0 && max_div < 1.0
          n_gaps = 0; n_opens = 0
          cg.scan(/(\d+)([MID])/) { |m| if m[2] != "M"
            n_gaps += m[1].to_i; n_opens += 1
          end }
          n_mm = t[10].to_i - t[9].to_i - n_gaps
          n_diff2 = n_mm + n_opens
          next if n_diff2.to_f / (n_diff2 + t[9].to_i) > max_div
        end

        strand = t[4]
        # Build sorted event list: each region contributes start and end events
        a = [] of {Int32, Int32, Int32, Int32} # {query_pos, type(0=st,1=en), reg_idx, ref_pos(-2=unset)}
        r = Array({Int32, Int32}).new(regs.size, {-2, -2})
        regs.each_with_index do |reg, ri|
          s = reg[0]; e = reg[1]
          if strand == "+"
            a << {s, 0, ri, -2}; a << {e - 1, 1, ri, -2}
          else
            a << {t[1].to_i - e, 0, ri, -2}; a << {t[1].to_i - s - 1, 1, ri, -2}
          end
        end
        a.sort_by! { |ev| ev[0] }

        # Walk CIGAR to map query positions to reference
        k = 0; rx = t[7].to_i
        y = strand == "+" ? t[2].to_i : t[1].to_i - t[3].to_i
        updated_a = a.map { |ev| [ev[0], ev[1], ev[2], ev[3]] }

        cg.scan(/(\d+)([MID])/) do |cm|
          len = cm[1].to_i; op = cm[2][0]
          if op == 'D'
            rx += len; next
          end
          while k < updated_a.size && updated_a[k][0] < y
            k += 1
          end
          j = k
          while j < updated_a.size
            break if updated_a[j][0] >= y + len
            if y <= updated_a[j][0] && updated_a[j][0] < y + len
              updated_a[j][3] = op == 'M' ? rx + (updated_a[j][0] - y) : rx
            end
            j += 1
          end
          y += len; rx += len if op == 'M'
        end

        updated_a.each do |ev|
          ri = ev[2]
          if ev[1] == 0
            r[ri] = {ev[3], r[ri][1]}
          else
            r[ri] = {r[ri][0], ev[3] + 1}
          end
        end
        r.each_with_index do |coords, ri|
          name = "#{t[0]}_#{regs[ri][0]}_#{regs[ri][1]}"
          r0 = coords[0]; r1 = coords[1]
          name += "_t5" if r0 < 0; r0 = t[7].to_i if r0 < 0
          name += "_t3" if r1 < 0; r1 = t[8].to_i if r1 < 0
          puts [t[5], r0, r1, name, 0, strand].join('\t')
        end
      end
    end
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # call — call variants from sorted PAF with cs tag
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_call(args : Array(String)) : Int32
    min_cov_len = 10000; min_var_len = 50000; gap_thres = 50; gap_thres_long = 1000
    min_mapq = 5; sample_name = "sample"
    fa_fn : String? = nil
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-l"; i += 1; min_cov_len = args[i].to_i
      when "-L"; i += 1; min_var_len = args[i].to_i
      when "-g"; i += 1; gap_thres = args[i].to_i
      when "-q"; i += 1; min_mapq = args[i].to_i
      when "-f"; i += 1; fa_fn = args[i]
      when "-s"; i += 1; sample_name = args[i]
      when "-h", "--help"
        STDERR.puts "Usage: sort -k6,6 -k8,8n <with-cs.paf> | paftools call [options] -"
        return 0
      else rest << args[i]
      end
      i += 1
    end
    if rest.empty?
      STDERR.puts "Usage: sort -k6,6 -k8,8n <with-cs.paf> | paftools call [options] -"; return 1
    end

    fa : Hash(String, String)? = nil
    fa_lens : Array({String, Int32})? = nil
    is_vcf = false
    if (fn = fa_fn)
      fa_data, lens_data = read_fasta(fn)
      fa = fa_data; fa_lens = lens_data; is_vcf = true
    end

    tot_len = 0_i64
    n_sub = [0_i64, 0_i64, 0_i64] # total, ts, tv
    n_ins = Array.new(5, 0_i64); n_del = Array.new(5, 0_i64)

    if is_vcf
      puts "##fileformat=VCFv4.1"
      fa_lens.try &.each { |name, l| puts "##contig=<ID=#{name},length=#{l}>" }
      puts "##INFO=<ID=QNAME,Number=1,Type=String,Description=\"Query name\">"
      puts "##INFO=<ID=QSTART,Number=1,Type=Integer,Description=\"Query start\">"
      puts "##INFO=<ID=QSTRAND,Number=1,Type=String,Description=\"Query strand\">"
      puts "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
      puts "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t#{sample_name}"
    end

    # out[i] = [ctg, x, x_end, cov, mapq, ref_allele, qry_allele, qname, qs, qe, strand]
    a_active = [] of Array(String) # active alignments: [ctg, start, end]
    out_buf = [] of Array(String)
    c1_ctg : String? = nil; c1_start = 0; c1_end = 0; c1_counted = false; c1_len = 0_i64

    flush_out = ->(upto_ctg : String?, upto_x : Int32) {
      while !out_buf.empty?
        o = out_buf[0]
        break if o[0] == upto_ctg && o[2].to_i > upto_x
        # count variant
        if o[3].to_i <= 1
          if o[5] == "-" && o[6] != "-" # insertion
            l = o[6].size
            bi = l == 1 ? 0 : l == 2 ? 1 : l < gap_thres ? 2 : l < gap_thres_long ? 3 : 4
            n_ins[bi] += 1
          elsif o[5] != "-" && o[6] == "-" # deletion
            l = o[5].size
            bi = l == 1 ? 0 : l == 2 ? 1 : l < gap_thres ? 2 : l < gap_thres_long ? 3 : 4
            n_del[bi] += 1
          elsif o[5] != "-" && o[6] != "-" # sub
            n_sub[0] += 1
            s = (o[5] + o[6]).downcase
            if s == "ag" || s == "ga" || s == "ct" || s == "tc"
              n_sub[1] += 1
            else
              n_sub[2] += 1
            end
          end
        end
        if is_vcf
          if o[3].to_i == 1 && !(o[5] == "-" && o[6] == "-")
            fasta = fa || Hash(String, String).new
            if o[5] != "-" && o[6] != "-" # SNP
              v = [o[0], (o[1].to_i + 1).to_s, ".", o[5].upcase, o[6].upcase, o[4], ".", "QNAME=#{o[7]};QSTART=#{o[8].to_i + 1};QSTRAND=#{o[10]}", "GT", "1/1"]
              puts v.join('\t')
            elsif o[1].to_i > 0
              ref_seq = fasta[o[0]]?
              if ref_seq
                ref_base = ref_seq[o[1].to_i - 1].upcase.to_s
                if o[5] == "-" # insertion
                  v = [o[0], o[1], ".", ref_base, ref_base + o[6].upcase, o[4], ".", "QNAME=#{o[7]};QSTART=#{o[8].to_i + 1};QSTRAND=#{o[10]}", "GT", "1/1"]
                else # deletion
                  v = [o[0], o[1], ".", ref_base + o[5].upcase, ref_base, o[4], ".", "QNAME=#{o[7]};QSTART=#{o[8].to_i + 1};QSTRAND=#{o[10]}", "GT", "1/1"]
                end
                puts v.join('\t')
              end
            end
          end
        else
          puts (["V"] + o).join('\t')
        end
        out_buf.shift
      end
    }

    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('\t', 12); next if t.size < 12 || t[5] == "*"
        (6..11).each { |j| t[j] = t[j].to_i.to_s }
        next if t[10].to_i < min_cov_len || t[11].to_i < min_mapq
        (1..3).each { |j| t[j] = t[j].to_i.to_s }

        ctg = t[5]; x = t[7].to_i; x_end = t[8].to_i
        query = t[0]; rev = t[4] == "-"; y = rev ? t[3].to_i : t[2].to_i
        cs : String? = nil; tp : String? = nil
        line.scan(/\t(\S\S:[AZif]):(\S+)/) do |tm|
          cs = tm[2] if tm[1] == "cs:Z"
          tp = tm[2] if tm[1] == "tp:A"
        end
        have_s1 = line.includes?("\ts1:i:"); have_s2 = line.includes?("\ts2:i:")
        next if have_s1 && !have_s2
        next if tp && (tp == "S" || tp == "i")

        # Update coverage-1 region tracking
        if ctg != c1_ctg || x >= c1_end
          if c1_counted && c1_end > c1_start
            c1_len += c1_end - c1_start
            puts (["R", c1_ctg, c1_start, c1_end].map(&.to_s)).join('\t') unless is_vcf
          end
          c1_ctg = ctg; c1_start = x; c1_end = x_end; c1_counted = t[10].to_i >= min_var_len
        elsif x_end > c1_end
          if c1_counted && x > c1_start
            c1_len += x - c1_start
            puts ["R", c1_ctg, c1_start, x].join('\t') unless is_vcf
          end
          c1_start = c1_end; c1_end = x_end; c1_counted = t[10].to_i >= min_var_len
        elsif x_end > c1_start
          if c1_counted && x > c1_start
            c1_len += x - c1_start
            puts ["R", c1_ctg, c1_start, x].join('\t') unless is_vcf
          end
          c1_start = x_end
        end

        flush_out.call(ctg, x)

        # Update coverage of buffered variants
        out_buf.each do |o|
          o[3] = (o[3].to_i + 1).to_s if o[1].to_i >= x && o[2].to_i <= x_end
        end

        # Drop stale active alignments
        a_active.reject! { |aa| aa[0] != ctg || aa[2].to_i <= x }

        # Core: parse cs tag for variants
        if t[10].to_i >= min_var_len && (cs_val = cs)
          tot_len += t[10].to_i
          rx = x; qy = y
          cs_val.scan(/([:=*+\-])([A-Za-z]+|\d+)/) do |cm|
            op = cm[1][0]; seq = cm[2]
            cov = 1
            if op == '*' || op == '+' || op == '-'
              a_active.each { |aa| cov += 1 if aa[2].to_i > rx }
            end
            case op
            when '=', ':'
              l = op == '=' ? seq.size : seq.to_i
              qy = rev ? qy - l : qy + l
              rx += l
            when '*'
              qs_pos = rev ? qy - 1 : qy; qe_pos = rev ? qy : qy + 1
              qy = rev ? qy - 1 : qy + 1
              br = seq[0].to_s; bq = seq[1].to_s
              unless br == "n" || bq == "n"
                out_buf << [ctg, rx.to_s, (rx + 1).to_s, cov.to_s, t[11], br, bq,
                            query, qs_pos.to_s, qe_pos.to_s, rev ? "-" : "+"]
              end
              rx += 1
            when '+'
              l = seq.size
              qs_pos = rev ? qy - l : qy; qe_pos = rev ? qy : qy + l
              qy = rev ? qy - l : qy + l
              out_buf << [ctg, rx.to_s, rx.to_s, cov.to_s, t[11], "-", seq,
                          query, qs_pos.to_s, qe_pos.to_s, rev ? "-" : "+"]
            when '-'
              l = seq.size
              out_buf << [ctg, rx.to_s, (rx + l).to_s, cov.to_s, t[11], seq, "-",
                          query, qy.to_s, qy.to_s, rev ? "-" : "+"]
              rx += l
            end
          end
        end
        a_active << [ctg, x.to_s, x_end.to_s]
      end
    end

    if c1_counted && c1_end > c1_start
      c1_len += c1_end - c1_start
      puts ["R", c1_ctg, c1_start, c1_end].join('\t') unless is_vcf
    end
    flush_out.call(nil, 0)

    STDERR.puts "#{c1_len} reference bases covered by exactly one contig"
    STDERR.puts "#{n_sub[0]} substitutions; ts/tv = #{n_sub[2] > 0 ? (n_sub[1].to_f/n_sub[2]).round(3) : "N/A"}"
    STDERR.puts "#{n_del[0]} 1bp deletions"; STDERR.puts "#{n_ins[0]} 1bp insertions"
    STDERR.puts "#{n_del[1]} 2bp deletions"; STDERR.puts "#{n_ins[1]} 2bp insertions"
    STDERR.puts "#{n_del[2]} [3,#{gap_thres}) deletions"; STDERR.puts "#{n_ins[2]} [3,#{gap_thres}) insertions"
    STDERR.puts "#{n_del[3]} [#{gap_thres},#{gap_thres_long}) deletions"; STDERR.puts "#{n_ins[3]} [#{gap_thres},#{gap_thres_long}) insertions"
    STDERR.puts "#{n_del[4]} >=#{gap_thres_long} deletions"; STDERR.puts "#{n_ins[4]} >=#{gap_thres_long} insertions"
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # gff2bed — convert GTF/GFF3 to BED12 (or junction BED with -j)
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_gff2bed(args : Array(String)) : Int32
    is_short = false; print_junc = false; output_gene = false; ens_canon_only = false
    fn_ucsc : String? = nil
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-s"          ; is_short = true
      when "-j"          ; print_junc = true
      when "-G"          ; output_gene = true
      when "-e"          ; ens_canon_only = true
      when "-u"          ; i += 1; fn_ucsc = args[i]
      when "-h", "--help"; STDERR.puts "Usage: paftools gff2bed [-sjGe] [-u fai] <in.gff>"; return 0
      else                 rest << args[i]
      end
      i += 1
    end
    if rest.empty?
      STDERR.puts "Usage: paftools gff2bed [-sjGe] [-u fai] <in.gff>"; return 1
    end

    ens2ucsc = Hash(String, String).new
    if (fai = fn_ucsc)
      File.open(fai) do |f|
        f.each_line(chomp: true) do |line|
          t = line.split('\t'); s = t[0]
          if s =~ /_(random|alt|decoy)$/
            s = s.gsub(/_(random|alt|decoy)$/, "").gsub(/^chr\S+_/, "")
          else
            s = s.gsub(/^chrUn_/, "")
          end
          s = s.gsub(/v(\d+)/, ".\\1")
          ens2ucsc[s] = t[0] if s != t[0]
        end
      end
    end

    colors = {"protein_coding" => "0,128,255", "mRNA" => "0,128,255",
              "lincRNA" => "0,192,0", "snRNA" => "0,192,0",
              "miRNA" => "0,192,0", "misc_RNA" => "0,192,0"}

    print_bed12 = ->(exons : Array(Array(String)), cds_st : Int32, cds_en : Int32) {
      return if exons.empty?
      name = is_short ? "#{exons[0][7]}|#{exons[0][5]}" : exons[0][4..6].join("|")
      sorted = exons.sort_by { |e| e[1].to_i }
      if print_junc
        (1...sorted.size).each do |j|
          puts [sorted[j][0], sorted[j - 1][2], sorted[j][1], name, 1000, sorted[j][3]].join('\t')
        end
        return
      end
      st = sorted[0][1].to_i; en = sorted.last[2].to_i
      cds_st = st if cds_st == (1 << 30); cds_en = en if cds_en == 0
      sizes = sorted.map { |e| e[2].to_i - e[1].to_i }
      starts = sorted.map { |e| e[1].to_i - st }
      color = colors[sorted[0][5]]? || "196,196,196"
      puts [sorted[0][0], st, en, name, 1000, sorted[0][3], cds_st, cds_en, color,
            sorted.size, sizes.join(",") + ",", starts.join(",") + ","].join('\t')
    }

    re_gtf = /\b(transcript_id|transcript_type|transcript_biotype|gene_name|gene_id|gbkey|transcript_name|tag) "([^"]+)";/
    re_gff3 = /\b(transcript_id|transcript_type|transcript_biotype|gene_name|gene_id|gbkey|transcript_name)=([^;]+)/
    re_gtf_gene = /\b(gene_id|gene_type|gene_name) "([^;]+)";/
    re_gff3_gene = /\b(gene_id|gene_type|source_gene|gene_biotype|gene_name)=([^;]+);/

    exons = [] of Array(String); cds_st = 1 << 30; cds_en = 0; last_id : String? = nil

    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('\t')
        next if t.size < 9 || t[0].starts_with?('#')

        if fn_ucsc
          t[0] = ens2ucsc[t[0]]? || t[0]
          puts t.join('\t'); next
        end

        if output_gene
          next if t[2] != "gene"
          id = nil; type = ""; name = "N/A"
          line.scan(re_gtf_gene) { |m| id = m[2] if m[1] == "gene_id"; type = m[2] if m[1] == "gene_type"; name = m[2] if m[1] == "gene_name" }
          line.scan(re_gff3_gene) { |m| id = m[2] if m[1] == "gene_id" || m[1] == "source_gene"; type = m[2] if m[1] == "gene_type" || m[1] == "gene_biotype"; name = m[2] if m[1] == "gene_name" }
          puts [t[0], t[3].to_i - 1, t[4], "#{id}|#{type}|#{name}", 1000, t[6]].join('\t')
          next
        end

        next if t[2] != "CDS" && t[2] != "exon"
        t[3] = (t[3].to_i - 1).to_s; t[4] = t[4].to_i.to_s

        id = nil; type = ""; name = "N/A"; biotype = ""
        tname = "N/A"; ens_canonical = false
        line.scan(re_gtf) { |m|
          case m[1]
          when "transcript_id"              ; id = m[2]
          when "transcript_type"            ; type = m[2]
          when "transcript_biotype", "gbkey"; biotype = m[2]
          when "gene_name", "gene_id"       ; name = m[2]
          when "transcript_name"            ; tname = m[2]
          when "tag"                        ; ens_canonical = true if m[2] == "Ensembl_canonical"
          end
        }
        line.scan(re_gff3) { |m|
          case m[1]
          when "transcript_id"              ; id = m[2]
          when "transcript_type"            ; type = m[2]
          when "transcript_biotype", "gbkey"; biotype = m[2]
          when "gene_name", "gene_id"       ; name = m[2]
          when "transcript_name"            ; tname = m[2]
          end
        }
        next if ens_canon_only && !ens_canonical
        type = biotype if type.empty? && !biotype.empty?
        next unless id
        if id != last_id
          print_bed12.call(exons, cds_st, cds_en)
          exons = [] of Array(String); cds_st = 1 << 30; cds_en = 0; last_id = id
        end
        chr = ens2ucsc[t[0]]? || t[0]
        chr = chr.gsub(/([A-Z]+\d+)\.(\d+)/, "chrUn_\\1v\\2") if fn_ucsc && chr =~ /^[A-Z]+\d+\.\d+$/
        if t[2] == "CDS"
          cds_st = [cds_st, t[3].to_i].min; cds_en = [cds_en, t[4].to_i].max
        else
          exons << [chr, t[3], t[4], t[6], id || "", type, name, tname]
        end
      end
    end
    print_bed12.call(exons, cds_st, cds_en) if last_id
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # gff2junc — extract splice junctions from GFF3
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_gff2junc(args : Array(String)) : Int32
    feat = "CDS"; rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-f"          ; i += 1; feat = args[i]
      when "-h", "--help"; STDERR.puts "Usage: paftools gff2junc [-f feature] <in.gff3>"; return 0
      else                 rest << args[i]
      end
      i += 1
    end
    if rest.empty?
      STDERR.puts "Usage: paftools gff2junc [-f feature] <in.gff3>"; return 1
    end

    process = ->(a : Array(Array(String))) {
      return if a.size < 2
      s = a.sort_by { |e| e[4].to_i }
      (1...s.size).each { |j| puts [s[j][1], s[j - 1][5], s[j][4], s[j][0], 0, s[j][7]].join('\t') }
    }

    buf = [] of Array(String)
    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('\t'); next if t.size < 9 || t[0][0] == '#'
        next if t[2].downcase != feat.downcase
        m = /\bParent=([^;]+)/.match(t[8])
        unless m
          STDERR.puts "Can't find Parent"; next
        end
        t[3] = (t[3].to_i - 1).to_s
        entry = [m[1], t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8]]
        if !buf.empty? && buf[0][0] != m[1]
          process.call(buf); buf.clear
        end
        buf << entry
      end
    end
    process.call(buf)
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # splice2bed — splice alignments (PAF/SAM) → BED12
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_splice2bed(args : Array(String)) : Int32
    keep_multi = false; fn_conv : String? = nil
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-m"          ; keep_multi = true
      when "-n"          ; i += 1; fn_conv = args[i]
      when "-h", "--help"; STDERR.puts "Usage: paftools splice2bed [-m] <in.paf>|<in.sam>"; return 0
      else                 rest << args[i]
      end
      i += 1
    end
    if rest.empty?
      STDERR.puts "Usage: paftools splice2bed [-m] <in.paf>|<in.sam>"; return 1
    end

    conv = if (fn = fn_conv)
      Hash(String, String).new.tap { |h|
        File.open(fn) { |f| f.each_line(chomp: true) { |l| t = l.split('\t'); h[t[0]] = t[1] } }
      }
    end

    colors = ["0,128,255", "255,0,0", "0,192,0"]
    print_lines = ->(a : Array(Array(String))) {
      return if a.empty?
      n_pri = a.count { |row| row[8] == "0" }
      a.each { |row| row[8] = "1" } if n_pri > 1
      STDERR.puts "Warning: #{a[0][3]} has no primary alignment" if n_pri == 0
      a.each do |row|
        next if !keep_multi && row[8] == "2"
        row[8] = colors[row[8].to_i]
        puts row.join('\t')
      end
    }

    buf = [] of Array(String)
    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        next if line.starts_with?('@')
        t = line.split('\t')
        qname = conv.try { |c| c[t[0]]? } || t[0]
        t[0] = qname

        is_pri = false; cigar : String? = nil; a1 : Array(String)? = nil

        if t.size >= 12 && (t[4] == "+" || t[4] == "-") # PAF
          t[12..].each do |tag|
            if tag.starts_with?("cg:Z:")
              cigar = tag[5..]
            elsif tag.starts_with?("s2:i:")
              is_pri = true
            end
          end
          a1 = [t[5], t[7], t[8], t[0], (t[9].to_f / t[10].to_f * 1000).to_i.to_s, t[4]]
        elsif t.size >= 10 # SAM
          flag = t[1].to_i
          next if (flag & 4) != 0
          cigar = t[5]; is_pri = (flag & 0x100) == 0
          if (flag & 1) != 0
            t[0] += "/#{(flag >> 6) & 3}"
          end
          te_val : String? = nil
          a1 = [t[2], (t[3].to_i - 1).to_s, te_val || "?", t[0], "1000", (flag & 16) != 0 ? "-" : "+"]
        else
          next
        end

        if !buf.empty? && buf[0][3] != t[0]
          print_lines.call(buf); buf.clear
        end

        next unless cigar
        x0 = 0; x = 0; bs = [] of Int32; bl = [] of Int32
        cigar.scan(/(\d+)([MIDNSHP=X])/) do |cm|
          l = cm[1].to_i; op = cm[2][0]
          case op
          when 'M', 'D'; x += l
          when 'N'
            bs << x0; bl << (x - x0); x0 = x + l; x += l
          when 'S', 'H' # ignore for BED
          end
        end
        bs << x0; bl << (x - x0)

        row = (a1 || next).dup
        if row[2] == "?"
          row[2] = (row[1].to_i + x).to_s
        end
        row.concat([row[1], row[2], is_pri ? "0" : "2",
                    bs.size.to_s, bl.map { |v| v.to_s }.join(",") + ",",
                    bs.map { |v| v.to_s }.join(",") + ","])
        buf << row
      end
    end
    print_lines.call(buf)
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # asmstat — assembly alignment statistics vs reference
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  private def self.n_stat(lens : Array(Int32), tot : Int64, quantile : Float64) : Int32
    sorted = lens.sort.reverse
    sum = 0_i64
    sorted.each do |l|
      return l if sum <= quantile * tot && sum + l > quantile * tot
      sum += l
    end
    0
  end

  private def self.aun_stat(lens : Array(Int32), tot : Int64) : Float64
    sorted = lens.sort.reverse
    x = 0_i64; y = 0.0
    sorted.each do |l|
      l2 = x + l <= tot ? l : (tot - x).to_i
      x += l; y += l2.to_f * l2.to_f / tot.to_f
      break if x >= tot
    end
    y
  end

  def self.cmd_asmstat(args : Array(String)) : Int32
    min_query_len = 0; min_seg_len = 10000; max_diff = 0.01
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-q"; i += 1; min_query_len = args[i].to_i
      when "-l"; i += 1; min_seg_len = args[i].to_i
      when "-d"; i += 1; max_diff = args[i].to_f
      when "-h", "--help"
        STDERR.puts "Usage: paftools asmstat [options] <ref.fa.fai> <asm1.paf> [...]"; return 0
      else rest << args[i]
      end
      i += 1
    end
    if rest.size < 2
      STDERR.puts "Usage: paftools asmstat [options] <ref.fa.fai> <asm1.paf> [...]"; return 1
    end

    ref_len = 0_i64
    File.open(rest[0]) { |f| f.each_line(chomp: true) { |l| ref_len += l.split('\t')[1].to_i64 } }

    labels = ["Length", "l_cov", "Rcov%", "Rdup%", "Qcov%", "NG75", "NG50", "NGA50", "AUNGA",
              "#breaks", "bp(#{min_seg_len},0)", "bp(#{min_seg_len},10k)"]
    rst = Array.new(labels.size) { [] of String }
    header = ["Metric"]
    n_asm = rest.size - 1

    n_asm.times do |ai|
      fn = rest[ai + 1]
      label = fn.sub(/\.paf(\.gz)?$/, "")
      header << label

      ref_blocks = [] of {String, Int32, Int32}
      qblock_len = [] of Int32
      bp = [] of {Int32, Int32} # {flank, gap}
      query = Hash(String, Int32).new
      n_breaks = 0; qcov = 0_i64
      last_qname : String? = nil
      qblocks = [] of {Int32, Int32, String, String, Int32, Int32}

      process_query = ->(qbs : Array({Int32, Int32, String, String, Int32, Int32})) {
        qbs_s = qbs.sort_by { |b| b[0] }
        st = -1; en = -1; q_cov = 0_i64; last_k : Int32? = nil; last_blen : Int32? = nil
        qbs_s.each_with_index do |b, k|
          blen = b[1] - b[0]
          if k > 0 && b[0] < qbs_s[k - 1][1]
            next if b[1] <= qbs_s[k - 1][1]
            blen = b[1] - qbs_s[k - 1][1]
          end
          qblock_len << blen
          if b[0] > en
            q_cov += en - st; st = b[0]; en = b[1]
          else
            en = [en, b[1]].max
          end
          if (lk = last_k) && (lb = last_blen)
            gap = 1_000_000_000
            if b[2] == qbs_s[lk][2] && b[3] == qbs_s[lk][3]
              g1 = b[0] - qbs_s[lk][1]
              g2 = b[3] == "+" ? b[4] - qbs_s[lk][5] : qbs_s[lk][4] - b[5]
              gap = (g1 - g2).abs
            end
            flank = [blen, lb].min
            bp << {flank, gap}
          end
          last_k = k; last_blen = blen
        end
        q_cov += en - st
        q_cov
      }

      open_in(fn) do |io|
        io.each_line(chomp: true) do |line|
          t = line.split('\t')
          next if t.size < 2
          ql = t[1].to_i
          next if ql < min_query_len
          query[t[0]] = ql
          next if t.size < 9 || t[5] == "*"
          next unless line =~ /\ttp:A:[PI]/
          cg = (m = /\tcg:Z:(\S+)/.match(line)) ? m[1] : nil
          nms = (m = /\tNM:i:(\d+)/.match(line)) ? m[1].to_i : nil
          if cg && nms
            n_m = 0; n_gapo = 0; n_gaps = 0
            cg.scan(/(\d+)([MID])/) { |cm| l = cm[1].to_i; if cm[2] == "M"
              n_m += l
            else
              n_gapo += 1; n_gaps += l
            end }
            raise "NM < gaps" if nms < n_gaps
            diff = (nms - n_gaps + n_gapo).to_f / (n_m + n_gapo)
            next if diff > max_diff
          end
          qs = t[2].to_i; qe = t[3].to_i; rs2 = t[7].to_i; re2 = t[8].to_i
          next if qe - qs < min_seg_len
          n_breaks += 1 if t[0] == last_qname
          if t[0] != last_qname
            if last_qname
              qcov += process_query.call(qblocks)
            end
            qblocks = [] of {Int32, Int32, String, String, Int32, Int32}
            last_qname = t[0]
          end
          ref_blocks << {t[5], rs2, re2}
          qblocks << {qs, qe, t[4], t[5], rs2, re2}
        end
      end
      qcov += process_query.call(qblocks) if last_qname

      asm_len = 0_i64; asm_lens = [] of Int32
      query.each { |_, l| asm_len += l; asm_lens << l }
      rst[0] << asm_len.to_s
      rst[5] << n_stat(asm_lens, ref_len, 0.75).to_s
      rst[6] << n_stat(asm_lens, ref_len, 0.50).to_s

      ref_blocks.sort_by! { |b| {b[0], b[1]} }
      l_cov = 0_i64
      last_ref : String? = nil; st2 = -1; en2 = -1
      ref_blocks.each do |b|
        if b[0] != last_ref || b[1] > en2
          l_cov += en2 - st2 if st2 >= 0 && last_ref
          last_ref = b[0]; st2 = b[1]; en2 = b[2]
        else
          en2 = [en2, b[2]].max
        end
      end
      l_cov += en2 - st2 if st2 >= 0
      rst[1] << l_cov.to_s
      rst[2] << "#{"%.2f" % (100.0 * l_cov / ref_len)}%"
      rst[4] << "#{"%.2f" % (100.0 * qcov / [asm_len, 1].max)}%"

      c1_ctg2 : String? = nil; c1_start2 = 0; c1_end2 = 0; c1_len2 = 0_i64
      ref_blocks.each do |b|
        if b[0] != c1_ctg2 || b[1] >= c1_end2
          c1_len2 += c1_end2 - c1_start2 if c1_end2 > c1_start2
          c1_ctg2 = b[0]; c1_start2 = b[1]; c1_end2 = b[2]
        elsif b[2] > c1_end2
          c1_len2 += b[1] - c1_start2 if b[1] > c1_start2
          c1_start2 = c1_end2; c1_end2 = b[2]
        elsif b[2] > c1_start2
          c1_len2 += b[1] - c1_start2 if b[1] > c1_start2
          c1_start2 = b[2]
        end
      end
      c1_len2 += c1_end2 - c1_start2 if c1_end2 > c1_start2
      rst[3] << (l_cov > 0 ? "#{"%.2f" % (100.0 * (l_cov - c1_len2) / l_cov)}%" : "0.00%")

      rst[7] << n_stat(qblock_len, ref_len, 0.50).to_s
      rst[8] << "%.0f" % aun_stat(qblock_len, ref_len)
      rst[9] << n_breaks.to_s
      count_bp = ->(min_blen : Int32, min_gap : Int32) { bp.count { |b| b[0] >= min_blen && b[1] >= min_gap } }
      rst[10] << count_bp.call(500, 0).to_s
      rst[11] << count_bp.call(500, 10000).to_s
    end

    puts header.join('\t')
    labels.size.times { |i| puts (["#{labels[i]}"] + rst[i]).join('\t') }
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # asmgene — gene completeness evaluation
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_asmgene(args : Array(String)) : Int32
    min_iden = 0.99_f64; min_cov = 0.99_f64; print_err = false; auto_only = false
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-i"; i += 1; min_iden = args[i].to_f
      when "-c"; i += 1; min_cov = args[i].to_f
      when "-e"; print_err = true
      when "-a"; auto_only = true
      when "-h", "--help"
        STDERR.puts "Usage: paftools asmgene [options] <ref-splice.paf> <asm-splice.paf> [...]"; return 0
      else rest << args[i]
      end
      i += 1
    end
    if rest.size < 2
      STDERR.puts "Usage: paftools asmgene [options] <ref-splice.paf> <asm-splice.paf> [...]"; return 1
    end

    process_qry = ->(a : Array(Array(Int32))) {
      b = a.select { |r| r[4] >= r[5] * min_iden }
      cnt = [0, 0.0_f64, 0]
      return cnt if b.empty?
      n_full = b.count { |r| r[3] - r[2] >= r[1] * min_cov }
      cnt[0] = n_full
      s = b.sort_by { |r| r[2] }
      l_cov2 = 0; st2 = s[0][2]; en2 = s[0][3]
      (1...s.size).each { |j| if s[j][2] <= en2
        en2 = [en2, s[j][3]].max
      else
        l_cov2 += en2 - st2; st2 = s[j][2]; en2 = s[j][3]
      end }
      l_cov2 += en2 - st2
      cnt[1] = b[0][1] > 0 ? l_cov2.to_f / b[0][1] : 0.0
      cnt[2] = b.size
      cnt
    }

    n_fn = rest.size
    gene = Hash(String, Array(Array(Int32 | Float64)?)).new
    header = [] of String
    refpos = Hash(String, Array(String)).new

    n_fn.times do |fi|
      fn = rest[fi]
      label = fn.sub(/\.paf(\.gz)?$/, "")
      header << label
      a = [] of Array(Int32)
      open_in(fn) do |io|
        io.each_line(chomp: true) do |line|
          t = line.split('\t')
          next if t.size < 12
          ql = t[1].to_i; qs2 = t[2].to_i; qe2 = t[3].to_i; mlen = t[9].to_i; blen = t[10].to_i
          refpos[t[0]] = [t[0], t[1], t[5], t[7], t[8]] if fi == 0
          gene[t[0]] ||= Array(Array(Int32 | Float64)?).new(n_fn, nil)
          if a.size > 0 && t[0] != a[0][0].to_s
            gene[a[0][0].to_s]?.try { |arr| arr[fi] = process_qry.call(a) }
            a = [] of Array(Int32)
          end
          a << [t[0].hash.to_i32, ql, qs2, qe2, mlen, blen]
        end
      end
      gene[a[0][0].to_s]?.try { |arr| arr[fi] = process_qry.call(a) } unless a.empty?
    end

    # Deduplicate reference: keep longest non-overlapping gene per locus
    gene_list = refpos.values.sort_by { |g| {g[2], g[3].to_i} }
    gene_nr = Hash(String, Bool).new; last_j = 0
    (1...gene_list.size).each do |j|
      if gene_list[j][2] != gene_list[last_j][2] || gene_list[j][3].to_i >= gene_list[last_j][4].to_i
        gene_nr[gene_list[last_j][0]] = true; last_j = j
      elsif gene_list[j][1].to_i > gene_list[last_j][1].to_i
        last_j = j
      end
    end
    gene_nr[gene_list[last_j][0]] = true unless gene_list.empty?

    col1 = ["full_sgl", "full_dup", "frag", "part50+", "part10+", "part10-", "dup_cnt", "dup_sum"]
    rst2 = Array.new(col1.size) { Array.new(n_fn, 0) }

    gene.each do |g, arr|
      next unless (ref_row = arr[0])
      ref_cnt = ref_row[0].as(Int32)
      next if ref_cnt != 1
      next unless gene_nr[g]
      next if auto_only && (rp = refpos[g]?) && rp[2] =~ /^(chr)?[XY]$/
      n_fn.times do |fi|
        if (ga = arr[fi]).nil?
          rst2[5][fi] += 1
          puts "M\t#{header[fi]}\t#{refpos[g]?.try(&.join("\t"))}" if print_err
        else
          cnt = ga
          case cnt[0].as(Int32)
          when 1 then rst2[0][fi] += 1
          when .>(1)
            rst2[1][fi] += 1
            puts "D\t#{header[fi]}\t#{refpos[g]?.try(&.join("\t"))}" if print_err
          else
            cov_frac = cnt[1].as(Float64)
            if cov_frac >= min_cov
              rst2[2][fi] += 1; puts "F\t#{header[fi]}\t#{refpos[g]?.try(&.join("\t"))}" if print_err
            elsif cov_frac >= 0.5
              rst2[3][fi] += 1
            elsif cov_frac >= 0.1
              rst2[4][fi] += 1
            else
              rst2[5][fi] += 1; puts "0\t#{header[fi]}\t#{refpos[g]?.try(&.join("\t"))}" if print_err
            end
          end
        end
      end
    end

    puts (["H", "Metric"] + header).join('\t')
    col1.size.times { |k| puts (["X", col1[k]] + rst2[k].map(&.to_s)).join('\t') }
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # bedcov — compute bases in BED regions covered by target BED
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_bedcov(args : Array(String)) : Int32
    print_len = false; to_merge = true
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-p"; print_len = true
      when "-d"; to_merge = false
      when "-h", "--help"
        STDERR.puts "Usage: paftools bedcov [-pd] <regions.bed> <target.bed>"; return 0
      else rest << args[i]
      end
      i += 1
    end
    if rest.size < 2
      STDERR.puts "Usage: paftools bedcov [-pd] <regions.bed> <target.bed>"; return 1
    end

    # Read regions BED and build interval index per chromosome
    reg_ivs = Hash(String, Array({Int32, Int32})).new
    reg_idx = Hash(String, Array(Int32)).new
    File.open(rest[0]) do |f|
      f.each_line(chomp: true) do |line|
        t = line.split('\t'); next if t.size < 3
        chr = t[0]; bst = t[1].to_i; ben = t[2].to_i
        reg_ivs[chr] ||= [] of {Int32, Int32}
        if t.size >= 12 && t[9] =~ /^\d+/
          n = t[9].to_i; sz = t[10].split(','); sts = t[11].split(',')
          n.times { |j| reg_ivs[chr] << {bst + sts[j].to_i, bst + sts[j].to_i + sz[j].to_i} }
        else
          reg_ivs[chr] << {bst, ben}
        end
      end
    end
    reg_ivs.each do |chr, ivs|
      intv_sort(ivs); intv_merge(ivs) if to_merge
      reg_idx[chr] = intv_build(ivs)
    end

    tot_len = 0_i64; hit_len = 0_i64

    open_in(rest[1]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('\t'); next if t.size < 3
        chr = t[0]; bst = t[1].to_i; ben = t[2].to_i
        segs = [] of {Int32, Int32}
        if t.size >= 12 && t[9] =~ /^\d+/
          n = t[9].to_i; sz = t[10].split(','); sts = t[11].split(',')
          n.times { |j| segs << {bst + sts[j].to_i, bst + sts[j].to_i + sz[j].to_i} }
        else
          segs << {bst, ben}
        end

        feat_len = segs.sum { |s| s[1] - s[0] }.to_i64
        tot_len += feat_len

        next unless reg_ivs.has_key?(chr)
        ivs = reg_ivs[chr]; idx = reg_idx[chr]
        overlap = [] of {Int32, Int32}
        segs.each do |seg|
          intv_ovlp(ivs, idx, seg[0], seg[1]).each do |ov|
            overlap << {[ov[0], seg[0]].max, [ov[1], seg[1]].min}
          end
        end
        feat_hit = overlap.empty? ? 0_i64 : cov_len(overlap)
        hit_len += feat_hit
        puts (["F", t[0..3].join('\t'), feat_len, feat_hit]).join('\t') if print_len
      end
    end

    STDERR.puts "# target bases: #{tot_len}"
    STDERR.puts "# target bases overlapping regions: #{hit_len} (#{"%.2f" % (tot_len > 0 ? 100.0*hit_len/tot_len : 0.0)}%)"
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # vcfpair — pair diploid VCF haplotypes
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_vcfpair(args : Array(String)) : Int32
    is_male = false; sample = "syndip"; hgver : String? = nil
    par = {"37" => [{0, 2699520}, {154931043, 155260560}]}
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-m"; is_male = true
      when "-s"; i += 1; sample = args[i]
      when "-g"; i += 1; hgver = args[i]
      when "-h", "--help"
        STDERR.puts "Usage: paftools vcfpair [-m] [-s sample] [-g 37] <in.pair.vcf>"; return 0
      else rest << args[i]
      end
      i += 1
    end
    if rest.empty?
      STDERR.puts "Usage: paftools vcfpair [-m] [-s sample] [-g 37] <in.pair.vcf>"; return 1
    end
    re_ctg = is_male ? /^(chr)?([0-9]+|X|Y)$/ : /^(chr)?([0-9]+|X)$/

    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        if line.starts_with?('#')
          next if line =~ /^##(source|reference)=/
          if (m = /^##contig=.*ID=([^\s,]+)/.match(line))
            next unless re_ctg.match(m[1])
          elsif line.starts_with?("#CHROM")
            t = line.split('\t'); t.delete_at(t.size - 1); t[-1] = sample
            puts "##FILTER=<ID=HET1,Description=\"Heterozygous in the first haplotype\">"
            puts "##FILTER=<ID=HET2,Description=\"Heterozygous in the second haplotype\">"
            puts "##FILTER=<ID=GAP1,Description=\"Uncalled in the first haplotype\">"
            puts "##FILTER=<ID=GAP2,Description=\"Uncalled in the second haplotype\">"
            line = t.join('\t')
          end
          puts line; next
        end
        t = line.split('\t')
        next unless re_ctg.match(t[0])
        gt : String? = nil; ad : Array(Int32)? = nil
        filters = [] of String; ht = [nil, nil] of String?
        [0, 1].each do |hi|
          m = /^(\.|[0-9]+)\/(\.|[0-9]+):(\S+)/.match(t[9 + hi])
          unless m
            STDERR.puts "Malformatted VCF"; return 1
          end
          s = m[3].split(',')
          the_ad = (ad ||= Array.new(s.size, 0))
          s.each_with_index { |v, j| the_ad[j] += v.to_i }
          if m[1] == '.'
            filters << "GAP#{hi + 1}"; ht[hi] = "."
          elsif m[1] != m[2]
            filters << "HET#{hi + 1}"; ht[hi] = "."
          else
            ht[hi] = m[1]
          end
        end
        t.delete_at(t.size - 1)
        hap = 0; st2 = t[1].to_i; en2 = st2 + t[3].size
        if is_male
          if t[0] =~ /^(chr)?X/
            if (hv = hgver) && (pars = par[hv]?)
              in_par = pars.any? { |r| r[0] <= st2 && en2 <= r[1] }
              hap = in_par ? 0 : 2
            end
          elsif t[0] =~ /^(chr)?Y/
            hap = 1
          end
        end
        if hap > 0 && filters.size == 1
          filters.clear if (hap == 2 && filters[0] == "GAP1") || (hap == 1 && filters[0] == "GAP2")
        end
        t[5] = "30"; t[6] = filters.empty? ? "." : filters.join(";")
        t << ht.map { |h| h || "." }.join("|") + ":" + (ad || [] of Int32).join(",")
        puts t.join('\t')
      end
    end
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # run — command dispatcher
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.run(argv : Array(String)) : Int32
    if argv.empty?
      STDERR.puts "Usage: paftools <command> [arguments]"
      STDERR.puts "Commands:"
      STDERR.puts "  view       convert PAF to BLAST-like (default), MAF, or LASTZ-cigar"
      STDERR.puts "  sam2paf    convert SAM to PAF"
      STDERR.puts "  delta2paf  convert MUMmer delta to PAF"
      STDERR.puts "  gff2bed    convert GFF/GTF to BED12"
      STDERR.puts "  gff2junc   convert GFF3 to junction BED"
      STDERR.puts "  splice2bed convert splice alignments (PAF/SAM) to BED12"
      STDERR.puts "  longcs2seq reconstruct sequences from long-cs PAF"
      STDERR.puts ""
      STDERR.puts "  stat       basic mapping statistics (PAF or SAM)"
      STDERR.puts "  asmstat    assembly alignment statistics"
      STDERR.puts "  asmgene    gene completeness evaluation"
      STDERR.puts "  liftover   lift BED coordinates via PAF CIGAR"
      STDERR.puts "  call       call variants from PAF with cs tag"
      STDERR.puts "  bedcov     compute coverage of BED regions"
      STDERR.puts "  vcfpair    pair diploid VCF haplotypes"
      STDERR.puts ""
      STDERR.puts "  version    print version"
      return 1
    end
    cmd = argv.shift
    case cmd
    when "stat"                 then cmd_stat(argv)
    when "view"                 then cmd_view(argv)
    when "sam2paf"              then cmd_sam2paf(argv)
    when "delta2paf"            then cmd_delta2paf(argv)
    when "liftover", "liftOver" then cmd_liftover(argv)
    when "call"                 then cmd_call(argv)
    when "longcs2seq"           then cmd_longcs2seq(argv)
    when "gff2bed"              then cmd_gff2bed(argv)
    when "gff2junc"             then cmd_gff2junc(argv)
    when "splice2bed"           then cmd_splice2bed(argv)
    when "asmstat"              then cmd_asmstat(argv)
    when "asmgene"              then cmd_asmgene(argv)
    when "bedcov"               then cmd_bedcov(argv)
    when "vcfpair"              then cmd_vcfpair(argv)
    when "version"              then puts VERSION; 0
    else                             STDERR.puts "Error: unknown command '#{cmd}'"; 1
    end
  end
end
