module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # ov-eval — estimate read-overlap sensitivity from a to-reference mapping
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_ov_eval(args : Array(String)) : Int32
    min_ovlp = 2000; min_frac = 0.95; min_mapq = 10
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-q"; i += 1; min_mapq = args[i].to_i
      when "-l"; i += 1; min_ovlp = args[i].to_i
      when "-f"; i += 1; min_frac = args[i].to_f
      when "-h", "--help"
        STDERR.puts "Usage: sort -k6,6 -k8,8n to-ref.paf | paftools ov-eval [options] - <ovlp.paf>"; return 0
      else rest << args[i]
      end
      i += 1
    end
    if rest.size < 2
      STDERR.puts "Usage: sort -k6,6 -k8,8n to-ref.paf | paftools ov-eval [options] - <ovlp.paf>"
      STDERR.puts "Options:"
      STDERR.puts "  -l INT     min overlap length [2000]"
      STDERR.puts "  -q INT     min mapping quality [10]"
      STDERR.puts "  -f FLOAT   min fraction of mapped length [0.95]"
      return 1
    end

    a = [] of {String, Int32, Int32, String}
    h = Hash(String, Int32).new

    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('\t')
        next if t[11].to_i < min_mapq
        is_pri = (12...t.size).any? { |j| t[j] == "tp:A:P" }
        next unless is_pri
        qlen = t[1].to_i; qs = t[2].to_i; qe = t[3].to_i
        ctg = t[5]; st = t[7].to_i; en = t[8].to_i
        next if qe - qs < min_ovlp || en - st < min_ovlp || (qe - qs).to_f64 / qlen < min_frac

        while a.size > 0
          break if a[0][0] == ctg && a[0][2] > st
          a.shift
        end
        a.each do |rec|
          next if rec[3] == t[0]
          len = [en, rec[2]].min - st
          if len >= min_ovlp
            key = rec[3] < t[0] ? "#{rec[3]}\t#{t[0]}" : "#{t[0]}\t#{rec[3]}"
            h[key] = len
          end
        end
        a << {ctg, st, en, t[0]}
      end
    end

    open_in(rest[1]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('\t')
        key = t[0] < t[5] ? "#{t[0]}\t#{t[5]}" : "#{t[5]}\t#{t[0]}"
        if (v = h[key]?) && v > 0
          h[key] = -v
        end
      end
    end

    n_ovlp = 0; n_missing = 0
    h.each_value do |v|
      n_ovlp += 1
      n_missing += 1 if v > 0
    end
    puts "#{n_ovlp} overlaps inferred from the reference mapping"
    puts "#{n_missing} missed by the read overlapper"
    rate = n_ovlp == 0 ? Float64::NAN : 100 * (1 - n_missing.to_f64 / n_ovlp)
    puts "#{js_fixed2(rate)}% sensitivity"
    0
  end
end
