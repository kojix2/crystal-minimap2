module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # asmstat — assembly alignment statistics vs reference
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  private def self.n_stat(lens : Array(Int32), tot : Int64, quantile : Float64) : Int32
    sorted = lens.sort.reverse!
    sum = 0_i64
    sorted.each do |len|
      return len if sum <= quantile * tot && sum + len > quantile * tot
      sum += len
    end
    0
  end

  private def self.aun_stat(lens : Array(Int32), tot : Int64) : Float64
    sorted = lens.sort.reverse!
    x = 0_i64; y = 0.0
    sorted.each do |len|
      l2 = x + len <= tot ? len : (tot - x).to_i
      x += len; y += l2.to_f * l2.to_f / tot.to_f
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
    File.open(rest[0]) { |file| file.each_line(chomp: true) { |line2| ref_len += line2.split('\t')[1].to_i64 } }

    labels = ["Length", "l_cov", "Rcov%", "Rdup%", "Qcov%", "NG75", "NG50", "NGA50", "AUNGA",
              "#breaks", "bp(#{min_seg_len},0)", "bp(#{min_seg_len},10k)"]
    rst = Array.new(labels.size) { [] of String }
    header = ["Metric"]
    n_asm = rest.size - 1

    n_asm.times do |asm_i|
      fn = rest[asm_i + 1]
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
            cg.scan(/(\d+)([MID])/) do |mat|
              l = mat[1].to_i
              if mat[2] == "M"
                n_m += l
              else
                n_gapo += 1; n_gaps += l
              end
            end
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
      query.each { |_, qlen| asm_len += qlen; asm_lens << qlen }
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
    labels.size.times { |k| puts (["#{labels[k]}"] + rst[k]).join('\t') }
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
      b = a.select { |row| row[4] >= row[5] * min_iden }
      cnt = [0, 0.0_f64, 0]
      return cnt if b.empty?
      n_full = b.count { |row| row[3] - row[2] >= row[1] * min_cov }
      cnt[0] = n_full
      s = b.sort_by { |row| row[2] }
      l_cov2 = 0; st2 = s[0][2]; en2 = s[0][3]
      (1...s.size).each do |j|
        if s[j][2] <= en2
          en2 = [en2, s[j][3]].max
        else
          l_cov2 += en2 - st2; st2 = s[j][2]; en2 = s[j][3]
        end
      end
      l_cov2 += en2 - st2
      cnt[1] = b[0][1] > 0 ? l_cov2.to_f / b[0][1] : 0.0
      cnt[2] = b.size
      cnt
    }

    n_fn = rest.size
    gene = Hash(String, Array(Array(Int32 | Float64)?)).new
    header = [] of String
    refpos = Hash(String, Array(String)).new

    n_fn.times do |fn_i|
      fn = rest[fn_i]
      label = fn.sub(/\.paf(\.gz)?$/, "")
      header << label
      a = [] of Array(Int32)
      open_in(fn) do |io|
        io.each_line(chomp: true) do |line|
          t = line.split('\t')
          next if t.size < 12
          ql = t[1].to_i; qs2 = t[2].to_i; qe2 = t[3].to_i; mlen = t[9].to_i; blen = t[10].to_i
          refpos[t[0]] = [t[0], t[1], t[5], t[7], t[8]] if fn_i == 0
          gene[t[0]] ||= Array(Array(Int32 | Float64)?).new(n_fn, nil)
          if a.size > 0 && t[0] != a[0][0].to_s
            gene[a[0][0].to_s]?.try { |arr| arr[fn_i] = process_qry.call(a) }
            a = [] of Array(Int32)
          end
          a << [t[0].hash.to_i32, ql, qs2, qe2, mlen, blen]
        end
      end
      gene[a[0][0].to_s]?.try { |arr| arr[fn_i] = process_qry.call(a) } unless a.empty?
    end

    # Deduplicate reference: keep longest non-overlapping gene per locus
    gene_list = refpos.values.sort_by! { |gen| {gen[2], gen[3].to_i} }
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

    gene.each do |gname, arr|
      next unless ref_row = arr[0]
      ref_cnt = ref_row[0].as(Int32)
      next if ref_cnt != 1
      next unless gene_nr[gname]
      next if auto_only && (rp = refpos[gname]?) && rp[2] =~ /^(chr)?[XY]$/
      n_fn.times do |fn_i|
        if (ga = arr[fn_i]).nil?
          rst2[5][fn_i] += 1
          puts "M\t#{header[fn_i]}\t#{refpos[gname]?.try(&.join("\t"))}" if print_err
        else
          cnt = ga
          case cnt[0].as(Int32)
          when 1 then rst2[0][fn_i] += 1
          when .>(1)
            rst2[1][fn_i] += 1
            puts "D\t#{header[fn_i]}\t#{refpos[gname]?.try(&.join("\t"))}" if print_err
          else
            cov_frac = cnt[1].as(Float64)
            if cov_frac >= min_cov
              rst2[2][fn_i] += 1; puts "F\t#{header[fn_i]}\t#{refpos[gname]?.try(&.join("\t"))}" if print_err
            elsif cov_frac >= 0.5
              rst2[3][fn_i] += 1
            elsif cov_frac >= 0.1
              rst2[4][fn_i] += 1
            else
              rst2[5][fn_i] += 1; puts "0\t#{header[fn_i]}\t#{refpos[gname]?.try(&.join("\t"))}" if print_err
            end
          end
        end
      end
    end

    puts (["H", "Metric"] + header).join('\t')
    col1.size.times { |k| puts (["X", col1[k]] + rst2[k].map(&.to_s)).join('\t') }
    0
  end
end
