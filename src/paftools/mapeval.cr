module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # mapeval — evaluate mapping accuracy against simulated-read truth encoded
  # in read names (pbsim single-end / mason2 paired-end conventions)
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  RE_PBSIM_SE  = /^(\S+)!(\S+)!(\d+)!(\d+)!([+-])$/
  RE_MASON2_PE = /^(\S+)!(\S+)!(\d+)_(\d+)!(\d+)_(\d+)!([+-])([+-])\/([12])$/

  # Truth interval: [chr, st, en, strand]
  private alias MapevalTruth = {String, Int32, Int32, String}
  # Candidate alignment: [tname, tstart, tend, strand, mapq, score]
  private alias MapevalAln = {String, Int32, Int32, String, Int32, Int32}

  private def self.mapeval_is_correct(s : MapevalTruth, b : MapevalAln, ovlp_ratio : Float64) : Bool
    return false if s[0] != b[0] || s[3] != b[3]
    o : Int32; l : Int32
    if s[1] < b[1]
      return false if s[2] <= b[1]
      o = [s[2], b[2]].min - b[1]
      l = [s[2], b[2]].max - s[1]
    else
      return false if b[2] <= s[1]
      o = [s[2], b[2]].min - s[1]
      l = [s[2], b[2]].max - b[1]
    end
    o.to_f64 / l > ovlp_ratio
  end

  private def self.mapeval_count_err(qname : String, a : Array(MapevalAln), tot : Array(Int32), err : Array(Int32),
                                      mode : Int32, err_out_q : Int32, ovlp_ratio : Float64, cap_short_mapq : Bool)
    return if a.empty?

    s : MapevalTruth
    if m = RE_PBSIM_SE.match(qname) # pbsim single-end reads
      s = {m[2], m[3].to_i, m[4].to_i, m[5]}
    elsif m = RE_MASON2_PE.match(qname) # mason2 paired-end reads
      s = m[9] == "1" ? {m[2], m[3].to_i, m[5].to_i, m[7]} : {m[2], m[4].to_i, m[6].to_i, m[8]}
    else
      raise "Failed to parse simulated read names '#{qname}'"
    end

    if mode == 0 || mode == 1 # longest only or first only
      max_i = 0
      if mode == 0 # longest only
        max = 0
        a.each_with_index { |b, i| if b[5] > max; max = b[5]; max_i = i; end }
      end
      mapq = a[max_i][4]
      tot[mapq] += 1
      unless mapeval_is_correct(s, a[max_i], ovlp_ratio)
        puts ["E", qname, a[max_i].to_a.join('\t')].join('\t') if mapq >= err_out_q
        err[mapq] += 1
      end
    elsif mode == 2 # all primary mode
      aa = a.dup
      max_err_mapq = -1; max_mapq_local = 0; max_err_i = -1
      if cap_short_mapq
        max = 0; max_q = 0
        aa.each { |b| if b[5] > max; max = b[5]; max_q = b[4]; end }
        aa = aa.map { |b| {b[0], b[1], b[2], b[3], (max_q < b[4] ? max_q : b[4]), b[5]} }
      end
      aa.each_with_index do |b, i|
        max_mapq_local = max_mapq_local > b[4] ? max_mapq_local : b[4]
        unless mapeval_is_correct(s, b, ovlp_ratio)
          if b[4] > max_err_mapq
            max_err_mapq = b[4]; max_err_i = i
          end
        end
      end
      if max_err_mapq >= 0
        tot[max_err_mapq] += 1
        err[max_err_mapq] += 1
        puts ["E", qname, aa[max_err_i].to_a.join('\t')].join('\t') if max_err_mapq >= err_out_q
      else
        tot[max_mapq_local] += 1
      end
    end
  end

  def self.cmd_mapeval(args : Array(String)) : Int32
    max_mapq = 60; mode = 0; err_out_q = 256; ovlp_ratio = 0.1; cap_short_mapq = false
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-Q"; i += 1; err_out_q = args[i].to_i
      when "-r"; i += 1; ovlp_ratio = args[i].to_f
      when "-m"; i += 1; mode = args[i].to_i
      when "-c"; cap_short_mapq = true
      when "-h", "--help"
        STDERR.puts "Usage: paftools mapeval [options] <in.paf>|<in.sam>"; return 0
      else rest << args[i]
      end
      i += 1
    end
    if rest.empty?
      STDERR.puts "Usage: paftools mapeval [options] <in.paf>|<in.sam>"
      STDERR.puts "Options:"
      STDERR.puts "  -r FLOAT   mapping correct if overlap_length/union_length>FLOAT [#{ovlp_ratio}]"
      STDERR.puts "  -Q INT     print wrong mappings with mapQ>=INT [don't print]"
      STDERR.puts "  -m INT     0: eval the longest aln only; 1: first aln only; 2: all primary aln [0]"
      return 1
    end

    tot = Array(Int32).new(max_mapq + 1, 0)
    err = Array(Int32).new(max_mapq + 1, 0)
    n_unmapped : Int32? = nil
    re_cigar = /(\d+)([MIDSHN=X])/

    last : String? = nil
    a = [] of MapevalAln

    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        next if line[0]? == '@'
        t = line.split('\t')
        if t[4] == "+" || t[4] == "-" # PAF
          if last != t[0]
            mapeval_count_err(last.not_nil!, a, tot, err, mode, err_out_q, ovlp_ratio, cap_short_mapq) if last
            a = [] of MapevalAln
            last = t[0]
          end
          # secondary alignment in minimap2 PAF: has s1:i but not s2:i
          next if /\ts1:i:\d+/.matches?(line) && !/\ts2:i:\d+/.matches?(line)
          mapq = t[11].to_i
          mapq = max_mapq if mapq > max_mapq
          a << {t[5], t[7].to_i, t[8].to_i, t[4], mapq, t[9].to_i}
        else # SAM
          flag = t[1].to_i
          read_no = (flag >> 6) & 0x3
          qname = t[0]
          unless /\/[12]$/.matches?(qname)
            qname = "#{t[0]}/#{read_no}" if read_no == 1 || read_no == 2
          end
          if last != qname
            mapeval_count_err(last.not_nil!, a, tot, err, mode, err_out_q, ovlp_ratio, cap_short_mapq) if last
            a = [] of MapevalAln
            last = qname
          end
          next if (flag & 0x100) != 0 # secondary alignment
          if (flag & 0x4) != 0 || t[2] == "*" # unmapped
            n_unmapped = (n_unmapped || 0) + 1
            next
          end
          mapq = t[4].to_i
          mapq = max_mapq if mapq > max_mapq
          pos = t[3].to_i - 1; pos_end = pos
          n_gap = 0; mlen = 0
          t[5].scan(re_cigar) do |m|
            len = m[1].to_i
            case m[2]
            when "M", "X", "="
              pos_end += len; mlen += len
            when "I"
              n_gap += len
            when "D"
              n_gap += len; pos_end += len
            end
          end
          score = pos_end - pos
          if m = /\tNM:i:(\d+)/.match(line)
            nm = m[1].to_i
            score = mlen - (nm - n_gap) if nm >= n_gap
          end
          a << {t[2], pos, pos_end, (flag & 16) != 0 ? "-" : "+", mapq, score}
        end
      end
    end
    mapeval_count_err(last.not_nil!, a, tot, err, mode, err_out_q, ovlp_ratio, cap_short_mapq) if last

    sum_tot = 0; sum_err = 0; q_out = -1; sum_tot2 = 0; sum_err2 = 0
    q = max_mapq
    while q >= 0
      unless tot[q] == 0
        if q_out < 0 || err[q] > 0
          if q_out >= 0
            rate = sum_err2.to_f64 / sum_tot2
            puts ["Q", q_out, sum_tot, sum_err, js_fixed(rate, 9), sum_tot2].join('\t')
          end
          sum_tot = 0; sum_err = 0; q_out = q
        end
        sum_tot += tot[q]; sum_err += err[q]
        sum_tot2 += tot[q]; sum_err2 += err[q]
      end
      q -= 1
    end
    rate = sum_tot2 == 0 ? Float64::NAN : sum_err2.to_f64 / sum_tot2
    puts ["Q", q_out, sum_tot, sum_err, js_fixed(rate, 9), sum_tot2].join('\t')
    puts ["U", n_unmapped].join('\t') if n_unmapped
    0
  end
end
