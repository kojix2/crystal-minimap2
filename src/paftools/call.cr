module Paftools
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
end
