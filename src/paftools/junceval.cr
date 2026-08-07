module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # junceval / exoneval — evaluate predicted splice junctions / exons in a
  # PAF/SAM/BED alignment against a reference GTF/GFF3 annotation
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  CHR_ONLY_RE = /^(chr)?([0-9]+|X|Y)$/

  # Reads a GTF/GFF3 file and groups lines whose feature column is in
  # `features` by transcript_id. Returns {chr => [{tid, ivs}, ...]}.
  # `features` accepts multiple names so callers matching either-case
  # "CDS"/"cds" (exoneval's `use_cds` mode) group same-tid rows together
  # in one pass, exactly mirroring the JS `t[2] != "cds" && t[2] != "CDS"`
  # single-pass filter (rather than reading each case separately and
  # merging afterwards, which would split a transcript's intervals across
  # two distinct tid entries for the rare file mixing both cases).
  private def self.read_gtf_transcripts(fn : String, features : Array(String)) : Hash(String, Array({String, Array({Int32, Int32})}))
    tr = Hash(String, {String, Array({Int32, Int32})}).new
    open_in(fn) do |io|
      io.each_line(chomp: true) do |line|
        next if line.starts_with?('#')
        t = line.split('\t')
        next unless features.includes?(t[2])
        st = t[3].to_i - 1
        en = t[4].to_i
        m = /transcript_id "(\S+)"/.match(t[8])
        next unless m
        tid = m[1]
        tr[tid] ||= {t[0], [] of {Int32, Int32}}
        tr[tid][1] << {st, en}
      end
    end
    result = Hash(String, Array({String, Array({Int32, Int32})})).new
    tr.each do |tid, (chr, ivs)|
      result[chr] ||= [] of {String, Array({Int32, Int32})}
      result[chr] << {tid, ivs}
    end
    result
  end

  private def self.junceval_usage
    STDERR.puts "Usage: paftools junceval [options] <gene.gtf> <aln.sam>"
    STDERR.puts "Options:"
    STDERR.puts "  -l INT    tolerance of junction positions (0 for exact) [0]"
    STDERR.puts "  -p        print overlapping introns"
    STDERR.puts "  -e        print erroreous overlapping introns"
    STDERR.puts "  -c        only consider alignments to /^(chr)?([0-9]+|X|Y)$/"
    STDERR.puts "  -a        miniprot PAF as input"
    STDERR.puts "  -b        BED as input"
    STDERR.puts "  -1        only process the first alignment of each query"
  end

  def self.cmd_junceval(args : Array(String)) : Int32
    l_fuzzy = 0; print_ovlp = false; print_err_only = false; first_only = false
    chr_only = false; aa = false; is_bed = false
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-l"; i += 1; l_fuzzy = args[i].to_i
      when "-e"; print_err_only = true; print_ovlp = true
      when "-p"; print_ovlp = true
      when "-c"; chr_only = true
      when "-a"; aa = true
      when "-b"; is_bed = true
      when "-1"; first_only = true
      when "-h", "--help"; junceval_usage; return 0
      else rest << args[i]
      end
      i += 1
    end
    if rest.empty?
      junceval_usage
      return 1
    end

    tr = read_gtf_transcripts(rest[0], ["exon"])
    anno = Hash(String, Array({Int32, Int32})).new
    tr.each do |chr, transcripts|
      anno[chr] ||= [] of {Int32, Int32}
      transcripts.each do |(tid, ivs)|
        s = ivs.dup
        intv_sort(s)
        (0...(s.size - 1)).each do |j|
          if s[j][1] >= s[j + 1][0]
            STDERR.puts "WARNING: incorrect annotation for transcript #{tid} (#{s[j][1]} >= #{s[j + 1][0]})"
          end
          anno[chr] << {s[j][1], s[j + 1][0]}
        end
      end
    end
    anno_idx = Hash(String, Array(Int32)).new
    anno.each do |chr, e|
      next if e.empty?
      intv_sort(e)
      k = 0
      (1...e.size).each do |idx|
        if e[idx] != e[k]
          k += 1
          e[k] = e[idx]
        end
      end
      e.truncate(0, k + 1)
      anno_idx[chr] = intv_build(e)
    end

    n_pri = 0; n_unmapped = 0; n_mapped = 0
    n_sgl = 0; n_splice = 0; n_splice_hit = 0; n_splice_novel = 0

    aln_fn = rest.size < 2 || rest[1] == "-" ? "-" : rest[1]
    last_qname : String? = nil
    re_cigar = /(\d+)([MIDNSHP=XFGUV])/

    open_in(aln_fn) do |io|
      io.each_line(chomp: true) do |line|
        next if line.starts_with?('@')
        t = line.split('\t')
        t.shift if t[0] == "##PAF"
        qname = t[0]
        ctg_name : String; pos : Int32; cigar : String? = nil

        if is_bed
          ctg_name = t[0]; pos = t[1].to_i
          # JS original: `ctg_name = t[0], pos = parseInt(t[1]), cigar == null;`
          # — the trailing `cigar == null` is a no-op equality test (should
          # have been `cigar = null`), but cigar already defaults to null at
          # this point, so the typo has no observable effect; kept as nil here.
        elsif t.size > 4 && (t[4] == "+" || t[4] == "-" || t[4] == "*") # PAF
          ctg_name = t[5]; pos = t[7].to_i
          type = "P"
          (12...t.size).each do |j|
            if m = /^(tp:A|cg:Z):(\S+)/.match(t[j])
              if m[1] == "tp:A"
                type = m[2]
              else
                cigar = m[2]
              end
            end
          end
          next if type == "S" # secondary
        else # SAM
          ctg_name = t[2]; pos = t[3].to_i - 1; cigar = t[5]
          flag = t[1].to_i
          if (flag & 1) != 0
            if (flag & 0x40) != 0
              qname += "/1"
            elsif (flag & 0x80) != 0
              qname += "/2"
            end
          end
          next if (flag & 0x100) != 0 # secondary
        end

        next if chr_only && !CHR_ONLY_RE.matches?(ctg_name)
        next if first_only && last_qname == qname
        if ctg_name == "*"
          n_unmapped += 1
          next
        else
          n_pri += 1
          if last_qname != qname
            n_mapped += 1
            last_qname = qname
          end
        end

        intron = [] of {Int32, Int32}
        if is_bed
          intron << {pos, t[2].to_i}
        elsif aa
          tmp_junc = [] of {Int32, Int32}; tmp = 0
          cigar.try &.scan(re_cigar) do |m|
            len = m[1].to_i; op = m[2]
            case op
            when "N"
              tmp_junc << {tmp, tmp + len}; tmp += len
            when "U"
              tmp_junc << {tmp + 1, tmp + len - 2}; tmp += len
            when "V"
              tmp_junc << {tmp + 2, tmp + len - 1}; tmp += len
            when "M", "X", "=", "D"
              tmp += len * 3
            when "F", "G"
              tmp += len
            end
          end
          if t[4] == "+"
            tmp_junc.each { |j| intron << {pos + j[0], pos + j[1]} }
          elsif t[4] == "-"
            glen = t[8].to_i - t[7].to_i
            tmp_junc.reverse_each { |j| intron << {pos + (glen - j[1]), pos + (glen - j[0])} }
          end
        else
          cigar.try &.scan(re_cigar) do |m|
            len = m[1].to_i; op = m[2]
            case op
            when "N"
              intron << {pos, pos + len}; pos += len
            when "M", "X", "=", "D"
              pos += len
            end
          end
        end

        if intron.empty?
          n_sgl += 1
          next
        end
        n_splice += intron.size

        chr = anno[ctg_name]?
        if chr
          idx = anno_idx[ctg_name]
          intron.each_with_index do |iv, k|
            o = intv_ovlp(chr, idx, iv[0], iv[1])
            if !o.empty?
              hit = false
              o.each do |ov|
                st_diff = (iv[0] - ov[0]).abs
                en_diff = (iv[1] - ov[1]).abs
                if st_diff <= l_fuzzy && en_diff <= l_fuzzy
                  n_splice_hit += 1; hit = true
                end
                break if hit
              end
              if print_ovlp
                type = hit ? "C" : "P"
                next if hit && print_err_only
                x = "[" + o.map { |ov| "(#{ov[0]},#{ov[1]})" }.join(", ") + "]"
                puts [type, qname, k + 1, ctg_name, iv[0], iv[1], x].join('\t')
              end
            else
              n_splice_novel += 1
              puts ["N", qname, k + 1, ctg_name, iv[0], iv[1]].join('\t') if print_ovlp
            end
          end
        else
          n_splice_novel += intron.size
        end
      end
    end

    unless print_ovlp
      puts "# unmapped reads: #{n_unmapped}"
      puts "# mapped reads: #{n_mapped}"
      puts "# primary alignments: #{n_pri}"
      puts "# singletons: #{n_sgl}"
      puts "# predicted introns: #{n_splice}"
      puts "# non-overlapping introns: #{n_splice_novel}"
      # js_fixed matches JS toFixed's "NaN" output when n_splice==0 (Crystal
      # sprintf would print "nan" instead).
      puts "# correct introns: #{n_splice_hit} (#{js_fixed(n_splice_hit.to_f / n_splice * 100, 2)}%)"
    end
    0
  end

  private def self.exoneval_usage
    STDERR.puts "Usage: paftools exoneval [options] <gene.gtf> <aln.sam>"
    STDERR.puts "Options:"
    STDERR.puts "  -l INT    tolerance of junction positions (0 for exact) [0]"
    STDERR.puts "  -d        evaluate coding regions only (exon regions by default)"
    STDERR.puts "  -a        miniprot PAF as input (force -d)"
    STDERR.puts "  -p        print overlapping exons"
    STDERR.puts "  -e        print erroreous overlapping exons"
    STDERR.puts "  -c        only consider alignments to /^(chr)?([0-9]+|X|Y)$/"
    STDERR.puts "  -1        only process the first alignment of each query"
    STDERR.puts "  -f        skip the first exon in the miniprot mode"
    STDERR.puts "  -t        skip the first and the last exons"
    STDERR.puts "  -b        BED as input"
    STDERR.puts "  -s        compute base Sn and Sp (more memory)"
  end

  private def self.exoneval_merge_and_index(ex : Hash(String, Array({Int32, Int32})))
    ex.each do |chr, e|
      s = e.dup
      intv_sort(s)
      a = [] of {Int32, Int32}
      st = s[0][0]; en = s[0][1]
      (1...s.size).each do |i|
        if s[i][0] > en
          a << {st, en}
          st = s[i][0]; en = s[i][1]
        else
          en = [en, s[i][1]].max
        end
      end
      a << {st, en}
      ex[chr] = a
    end
  end

  private def self.exoneval_cal_sn(a0 : Hash(String, Array({Int32, Int32})), a1 : Hash(String, Array({Int32, Int32})))
    tot = 0_i64; cov = 0_i64
    a1.each do |chr, e1|
      e1.each { |iv| tot += iv[1] - iv[0] }
      e0 = a0[chr]?
      next unless e0
      idx0 = intv_build(e0)
      e1.each do |iv|
        o = intv_ovlp(e0, idx0, iv[0], iv[1])
        o.each do |ov|
          st = [iv[0], ov[0]].max
          en = [iv[1], ov[1]].min
          cov += en - st
        end
      end
    end
    {tot, cov}
  end

  def self.cmd_exoneval(args : Array(String)) : Int32
    l_fuzzy = 0; print_ovlp = false; print_err_only = false; first_only = false
    chr_only = false; aa = false; is_bed = false; use_cds = false; eval_base = false
    skip_start = false; skip_last = false
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-l"; i += 1; l_fuzzy = args[i].to_i
      when "-e"; print_err_only = true; print_ovlp = true
      when "-p"; print_ovlp = true
      when "-c"; chr_only = true
      when "-a"; aa = true; use_cds = true
      when "-b"; is_bed = true
      when "-1"; first_only = true
      when "-d"; use_cds = true
      when "-s"; eval_base = true
      when "-f"; skip_start = true
      when "-t"; skip_last = true; skip_start = true
      when "-h", "--help"; exoneval_usage; return 0
      else rest << args[i]
      end
      i += 1
    end
    if rest.empty?
      exoneval_usage
      return 1
    end

    STDERR.puts "Reading reference GTF..."
    tr = read_gtf_transcripts(rest[0], use_cds ? ["CDS", "cds"] : ["exon"])

    anno = Hash(String, Array({Int32, Int32})).new
    tr.each do |chr, transcripts|
      anno[chr] ||= [] of {Int32, Int32}
      transcripts.each do |(_tid, ivs)|
        s = ivs.dup
        intv_sort(s)
        s.each { |iv| anno[chr] << iv }
      end
    end
    anno_idx = Hash(String, Array(Int32)).new
    anno.each do |chr, e|
      next if e.empty?
      intv_sort(e)
      k = 0
      (1...e.size).each do |idx|
        if e[idx] != e[k]
          k += 1
          e[k] = e[idx]
        end
      end
      e.truncate(0, k + 1)
      anno_idx[chr] = intv_build(e)
    end

    n_pri = 0; n_unmapped = 0; n_mapped = 0
    n_exon = 0; n_exon_hit = 0; n_exon_novel = 0
    qexon = Hash(String, Array({Int32, Int32})).new

    aln_fn = rest.size < 2 || rest[1] == "-" ? "-" : rest[1]
    last_qname : String? = nil
    re_cigar = /(\d+)([MIDNSHP=XFGUV])/

    STDERR.puts "Evaluating alignments..."
    open_in(aln_fn) do |io|
      io.each_line(chomp: true) do |line|
        next if line.starts_with?('@')
        t = line.split('\t')
        t.shift if t[0] == "##PAF"
        qname = t[0]
        ctg_name : String; pos : Int32; cigar : String? = nil

        if is_bed
          ctg_name = t[0]; pos = t[1].to_i
        elsif t.size > 4 && (t[4] == "+" || t[4] == "-" || t[4] == "*") # PAF
          ctg_name = t[5]; pos = t[7].to_i
          type = "P"
          (12...t.size).each do |j|
            if m = /^(tp:A|cg:Z):(\S+)/.match(t[j])
              if m[1] == "tp:A"
                type = m[2]
              else
                cigar = m[2]
              end
            end
          end
          next if type == "S" # secondary
        else # SAM
          ctg_name = t[2]; pos = t[3].to_i - 1; cigar = t[5]
          flag = t[1].to_i
          next if (flag & 0x100) != 0 # secondary
        end

        next if chr_only && !CHR_ONLY_RE.matches?(ctg_name)
        next if first_only && last_qname == qname
        if ctg_name == "*"
          n_unmapped += 1
          next
        else
          n_pri += 1
          if last_qname != qname
            n_mapped += 1
            last_qname = qname
          end
        end

        exon = [] of {Int32, Int32}
        if is_bed
          exon << {pos, t[2].to_i}
        elsif aa
          tmp_exon = [] of {Int32, Int32}; tmp = 0; tmp_st = 0
          cigar.try &.scan(re_cigar) do |m|
            len = m[1].to_i; op = m[2]
            case op
            when "N"
              tmp_exon << {tmp_st, tmp}; tmp_st = tmp + len; tmp += len
            when "U"
              tmp_exon << {tmp_st, tmp + 1}; tmp_st = tmp + len - 2; tmp += len
            when "V"
              tmp_exon << {tmp_st, tmp + 2}; tmp_st = tmp + len - 1; tmp += len
            when "M", "X", "=", "D"
              tmp += len * 3
            when "F", "G"
              tmp += len
            end
          end
          tmp_exon << {tmp_st, tmp}
          if t[4] == "+"
            tmp_exon.each { |e| exon << {pos + e[0], pos + e[1]} }
          elsif t[4] == "-" # coords on the query strand for protein-to-genome; flip them
            glen = t[8].to_i - t[7].to_i
            tmp_exon.reverse_each { |e| exon << {pos + (glen - e[1]), pos + (glen - e[0])} }
          end
          exon.shift if skip_start
          exon.pop if skip_last
        else
          tmp_st = pos
          cigar.try &.scan(re_cigar) do |m|
            len = m[1].to_i; op = m[2]
            case op
            when "N"
              exon << {tmp_st, pos}; tmp_st = pos + len; pos += len
            when "M", "X", "=", "D"
              pos += len
            end
          end
          exon << {tmp_st, pos}
        end
        n_exon += exon.size

        chr = anno[ctg_name]?
        if chr
          idx = anno_idx[ctg_name]
          exon.each_with_index do |iv, k|
            if eval_base
              qexon[ctg_name] ||= [] of {Int32, Int32}
              qexon[ctg_name] << iv
            end
            o = intv_ovlp(chr, idx, iv[0], iv[1])
            if !o.empty?
              hit = false
              o.each do |ov|
                st_diff = (iv[0] - ov[0]).abs
                en_diff = (iv[1] - ov[1]).abs
                if st_diff <= l_fuzzy && en_diff <= l_fuzzy
                  n_exon_hit += 1; hit = true
                end
                break if hit
              end
              if print_ovlp
                type = hit ? "C" : "P"
                next if hit && print_err_only
                x = "[" + o.map { |ov| "(#{ov[0]},#{ov[1]})" }.join(", ") + "]"
                puts [type, qname, k + 1, ctg_name, iv[0], iv[1], x].join('\t')
              end
            else
              n_exon_novel += 1
              puts ["N", qname, k + 1, ctg_name, iv[0], iv[1]].join('\t') if print_ovlp
            end
          end
        else
          n_exon_novel += exon.size
        end
      end
    end

    unless print_ovlp
      puts "# unmapped reads: #{n_unmapped}"
      puts "# mapped reads: #{n_mapped}"
      puts "# primary alignments: #{n_pri}"
      puts "# predicted exons: #{n_exon}"
      puts "# non-overlapping exons: #{n_exon_novel}"
      # js_fixed matches JS toFixed's "NaN" output when n_exon==0.
      puts "# correct exons: #{n_exon_hit} (#{js_fixed(n_exon_hit.to_f / n_exon * 100, 2)}%)"
    end

    if eval_base
      STDERR.puts "Computing base Sn and Sp..."
      exoneval_merge_and_index(qexon)
      exoneval_merge_and_index(anno)
      sn = exoneval_cal_sn(qexon, anno)
      sp = exoneval_cal_sn(anno, qexon)
      # js_fixed matches JS toFixed's "NaN" output when the denominator is 0.
      puts "Base Sn: #{sn[1]} / #{sn[0]} = #{js_fixed(sn[1].to_f / sn[0] * 100, 2)}%"
      puts "Base Sp: #{sp[1]} / #{sp[0]} = #{js_fixed(sp[1].to_f / sp[0] * 100, 2)}%"
    end
    0
  end
end
