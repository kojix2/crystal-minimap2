module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # misjoin — detect candidate misjoins/inversions/gaps in an assembly-to-
  # reference PAF (assumes PAF is grouped/sorted by query name, col 0)
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_misjoin(args : Array(String)) : Int32
    min_seg_len = 1_000_000; max_gap = 1_000_000
    fn_cen : String? = nil; show_long = false; show_err = false; cen_ratio = 0.5
    n_diff = [0, 0]; n_gap = [0, 0]; n_inv = [0, 0]; n_inv_end = [0, 0]

    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-l"          ; i += 1; min_seg_len = parse_num(args[i])
      when "-g"          ; i += 1; max_gap = parse_num(args[i])
      when "-c"          ; i += 1; fn_cen = args[i]
      when "-r"          ; i += 1; cen_ratio = args[i].to_f
      when "-p"          ; show_long = true
      when "-e"          ; show_err = true
      when "-h", "--help"; STDERR.puts "Usage: paftools misjoin [options] <in.paf>"; return 0
      else                 rest << args[i]
      end
      i += 1
    end
    if rest.empty?
      STDERR.puts "Usage: paftools misjoin [options] <in.paf>"
      STDERR.puts "Options:"
      STDERR.puts "  -c FILE   BED for centromeres []"
      STDERR.puts "  -r FLOAT  count a centromeric event if overlap ratio > FLOAT [#{cen_ratio}]"
      STDERR.puts "  -l NUM    min alignment block length [1m]"
      STDERR.puts "  -g NUM    max gap size [1m]"
      STDERR.puts "  -e        output misjoins not involving centromeres"
      STDERR.puts "  -p        output long alignment blocks for debugging"
      return 1
    end

    cen = Hash(String, Array({Int32, Int32})).new
    if fn = fn_cen
      open_in(fn) do |file|
        file.each_line(chomp: true) do |line|
          t = line.split('\t')
          cen[t[0]] ||= [] of {Int32, Int32}
          cen[t[0]] << {t[1].to_i, t[2].to_i}
        end
      end
    end

    test_cen = ->(chr : String, st : Int32, en : Int32) {
      b = cen[chr]?
      if b.nil?
        false
      else
        len = 0
        b.each do |seg|
          if seg[0] < en && seg[1] > st
            s = [seg[0], st].max
            e = [seg[1], en].min
            len += e - s
          end
        end
        len >= (en - st) * cen_ratio
      end
    }

    test_cen_point = ->(chr : String, x : Int32) {
      b = cen[chr]?
      b.nil? ? false : b.any? { |seg| x >= seg[0] && x < seg[1] }
    }

    if show_err || show_long
      puts "C\tJ  inter-chromosomal misjoin"
      puts "C\tj  inter-chromosomal misjoin with both breakpoints ending in centromeres"
      puts "C\tG  long gap on the reference genome"
      puts "C\tg  long gap on the reference genome with both breakpoints ending in centromeres"
      puts "C\tM  closed inversion"
      puts "C"
    end

    process = ->(rows : Array(Array(String))) {
      # rows[k]: 0=qname 1=qlen 2=qs 3=qe 4=strand 5=tname 6=tlen 7=ts 8=te 9=nmatch 10=blen 11=mapq
      kept = [] of Array(String)
      rows.each do |row|
        blen = row[10].to_i
        kept << row if blen >= min_seg_len
      end
      if kept.size > 1
      a = kept.sort_by { |row| row[2].to_i }
      # JS: print(a[i].join("\t")) prints the FULL row (all columns, incl.
      # any trailing PAF tags), unlike the error-report prints below which
      # deliberately slice(0, 12). Do not truncate here.
      if show_long
        a.each { |row| puts row.join('\t') }
      end
      i = 1
      while i < a.size
        ov0 = test_cen.call(a[i - 1][5], a[i - 1][7].to_i, a[i - 1][8].to_i)
        ov1 = test_cen.call(a[i][5], a[i][7].to_i, a[i][8].to_i)
        end0 = test_cen_point.call(a[i - 1][5], a[i - 1][4] == "+" ? a[i - 1][8].to_i : a[i - 1][7].to_i)
        end1 = test_cen_point.call(a[i][5], a[i][4] == "+" ? a[i][7].to_i : a[i][8].to_i)

        if a[i - 1][5] != a[i][5] # different chr
          if ov0 || ov1
            n_diff[1] += 1
          elsif show_err
            label = (end0 && end1) ? "j" : "J"
            puts [label, a[i - 1][0..11].join('\t')].join('\t')
            puts [label, a[i][0..11].join('\t')].join('\t')
          end
          n_diff[0] += 1
        elsif a[i - 1][4] == a[i][4] # a gap
          dq = a[i][2].to_i - a[i - 1][3].to_i
          dr = a[i][4] == "+" ? a[i][7].to_i - a[i - 1][8].to_i : a[i - 1][7].to_i - a[i][8].to_i
          gap = (dr - dq).abs
          if gap > max_gap
            if ov0 || ov1
              n_gap[1] += 1
            elsif show_err
              label = (end0 && end1) ? "g" : "G"
              puts [label, a[i - 1][0..11].join('\t')].join('\t')
              puts [label, a[i][0..11].join('\t')].join('\t')
            end
            n_gap[0] += 1
          end
        elsif i + 1 < a.size && a[i + 1][4] == a[i - 1][4] # bracketed inversion
          if ov0 || ov1
            n_inv[1] += 1
          elsif show_err
            puts ["M", a[i - 1][0..11].join('\t')].join('\t')
            puts ["M", a[i][0..11].join('\t')].join('\t')
            puts ["M", a[i + 1][0..11].join('\t')].join('\t')
          end
          n_inv[0] += 1
          i += 1
        else # hanging inversion
          if ov0 || ov1
            n_inv_end[1] += 1
          end
          n_inv_end[0] += 1
        end
        i += 1
      end
      end
    }

    a = [] of Array(String)
    open_in(rest[0] == "-" ? "-" : rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('\t')
        if a.size > 0 && a[0][0] != t[0]
          process.call(a)
          a = [] of Array(String)
        end
        a << t
      end
    end
    process.call(a) if a.size > 0

    puts "# inter-chromosomal misjoins: #{n_diff.join(",")}"
    puts "# intra-chromosomal gaps: #{n_gap.join(",")}"
    puts "# candidate inversions in the middle: #{n_inv.join(",")}"
    puts "# candidate inversions at contig ends: #{n_inv_end.join(",")}"
    0
  end
end
