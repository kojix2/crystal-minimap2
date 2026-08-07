module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # vcfstat — tabulate substitution/indel counts from a VCF
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  TRANSITIONS = {"AG" => true, "GA" => true, "CT" => true, "TC" => true}

  def self.cmd_vcfstat(args : Array(String)) : Int32
    rest = [] of String
    args.each do |a|
      if a == "-h" || a == "--help"
        STDERR.puts "Usage: paftools vcfstat <in.vcf>"; return 0
      else
        rest << a
      end
    end

    sub = 0_i64; ts = 0_i64; tv = 0_i64
    ins = 0_i64; del = 0_i64
    ins1 = 0_i64; del1 = 0_i64; ins2 = 0_i64; del2 = 0_i64
    ins50 = 0_i64; del50 = 0_i64; ins1k = 0_i64; del1k = 0_i64
    ins7k = 0_i64; del7k = 0_i64; insinf = 0_i64; delinf = 0_i64

    src = rest.empty? ? "-" : rest[0]
    open_in(src) do |io|
      io.each_line(chomp: true) do |line|
        next if line.starts_with?('#')
        t = line.split('\t')
        alts = t[4].split(',')
        ref = t[3]
        alts.each do |a|
          # The JS original tests `a[0]=='<' || a[1]=='>'` to skip symbolic
          # ALT alleles (e.g. "<DEL>"). The second half of that OR looks like
          # a typo — checking the *second* character for '>' only matches a
          # 2-character allele string, never the closing bracket of a real
          # symbolic allele. The rest of this file (_paf_get_alen, ported as
          # Paftools.paf_get_alen) instead uses /^<\S+>$/ to recognize
          # symbolic alleles, so we use that same, presumably-intended check
          # here for consistency instead of replicating the typo.
          next if /^<\S+>$/.matches?(a)
          l = [ref.size, a.size].min
          (0...l).each do |j|
            if ref[j] != a[j]
              sub += 1
              if TRANSITIONS["#{ref[j]}#{a[j]}"]?
                ts += 1
              else
                tv += 1
              end
            end
          end
          d = a.size - ref.size
          if d > 0
            ins += 1
            if d == 1
              ins1 += 1
            elsif d == 2
              ins2 += 1
            elsif d < 50
              ins50 += 1
            elsif d < 1000
              ins1k += 1
            elsif d < 7000
              ins7k += 1
            else
              insinf += 1
            end
          elsif d < 0
            d = -d
            del += 1
            if d == 1
              del1 += 1
            elsif d == 2
              del2 += 1
            elsif d < 50
              del50 += 1
            elsif d < 1000
              del1k += 1
            elsif d < 7000
              del7k += 1
            else
              delinf += 1
            end
          end
        end
      end
    end

    puts "# substitutions: #{sub}"
    # js_fixed handles the tv==0 edge case the same way JS's toFixed does
    # (prints "NaN"/"Infinity" rather than Crystal sprintf's "nan"/"inf").
    puts "ts/tv: #{js_fixed(ts.to_f / tv, 3)}"
    puts "# insertions: #{ins}"
    puts "# 1bp insertions: #{ins1}"
    puts "# 2bp insertions: #{ins2}"
    puts "# [3,50) insertions: #{ins50}"
    puts "# [50,1000) insertions: #{ins1k}"
    puts "# [1000,7000) insertions: #{ins7k}"
    puts "# >=7000 insertions: #{insinf}"
    puts "# deletions: #{del}"
    puts "# 1bp deletions: #{del1}"
    puts "# 2bp deletions: #{del2}"
    puts "# [3,50) deletions: #{del50}"
    puts "# [50,1000) deletions: #{del1k}"
    puts "# [1000,7000) deletions: #{del7k}"
    puts "# >=7000 deletions: #{delinf}"
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # vcfsel — filter VCF records by SV length
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_vcfsel(args : Array(String)) : Int32
    min_l = 0; max_l = 1 << 30
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-l"          ; i += 1; min_l = args[i].to_i
      when "-L"          ; i += 1; max_l = args[i].to_i
      when "-h", "--help"; STDERR.puts "Usage: paftools vcfsel [options] <in.vcf>"; return 0
      else                 rest << args[i]
      end
      i += 1
    end
    if rest.empty?
      STDERR.puts "Usage: paftools vcfsel [options] <in.vcf>"
      return 1
    end

    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        if line.starts_with?('#')
          puts line; next
        end
        t = line.split('\t')
        st = t[1].to_i; en = st + t[3].size - 1
        if m = /(^|;)END=(\d+)/.match(t[7])
          en = m[2].to_i
        end
        if en < st
          STDERR.puts "END is smaller than POS: #{en} < #{st}"
          en = st
        end
        _alen, min_abs_diff, max_abs_diff = paf_get_alen(t)
        next if max_abs_diff < min_l || min_abs_diff > max_l
        puts line
      end
    end
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # sveval — evaluate a called SV VCF against a base/truth VCF
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  # {st, en, svlen, abslen}
  private alias SvRec = {Int32, Int32, Int32?, Int32}

  private def self.sveval_read_bed(fn : String) : {Hash(String, Array({Int32, Int32})), Hash(String, Array(Int32))}
    ivs = Hash(String, Array({Int32, Int32})).new
    File.open(fn) do |file|
      file.each_line(chomp: true) do |line|
        t = line.split('\t')
        ivs[t[0]] ||= [] of {Int32, Int32}
        ivs[t[0]] << {t[1].to_i, t[2].to_i}
      end
    end
    idx = Hash(String, Array(Int32)).new
    ivs.each do |chr, a|
      intv_sort(a); intv_merge(a)
      idx[chr] = intv_build(a)
    end
    {ivs, idx}
  end

  private def self.sveval_read_vcf(fn : String, bed : Hash(String, Array({Int32, Int32}))?,
                                    bed_idx : Hash(String, Array(Int32))?,
                                    min_flt : Int32, max_size : Int32) : Hash(String, Array(SvRec))
    v = Hash(String, Array(SvRec)).new
    open_in(fn) do |io|
      io.each_line(chomp: true) do |line|
        next if line.starts_with?('#')
        t = line.split('\t')
        next if bed && !bed.has_key?(t[0])
        next if t[4] == "<INV>" || t[4] == "<INVDUP>" # no inversion
        next if /[\[\]]/.matches?(t[4])                # no break points
        next if t[6] != "." && t[6] != "PASS"
        st = t[1].to_i - 1
        en = st + t[3].size
        alen, _min_abs, _max_abs = paf_get_alen(t)
        svlen = alen
        abslen = svlen ? (svlen > 0 ? svlen : -svlen) : 0
        next if abslen < min_flt || abslen > max_size
        if m = /(^|;)END=(\d+)/.match(t[7])
          en = m[2].to_i
        elsif svlen && svlen < 0
          en = st + (-svlen)
        end
        en = st if en < st
        if st == en
          st -= 1; en += 1
        end
        if bed && bed_idx
          b = bed[t[0]]?; bi = bed_idx[t[0]]?
          next if b.nil? || bi.nil? || intv_ovlp(b, bi, st, en).empty?
        end
        v[t[0]] ||= [] of SvRec
        v[t[0]] << {st, en, svlen, abslen}
      end
    end
    v.each_value { |a| a.sort_by! { |r| {r[0], r[1]} } }
    v
  end

  private def self.sveval_idx(a : Array(SvRec)) : Array(Int32)
    ivs = a.map { |r| {r[0], r[1]} }
    intv_build(ivs)
  end

  private def self.sveval_ovlp(a : Array(SvRec), idx : Array(Int32), st : Int32, en : Int32) : Array(SvRec)
    ivs = a.map { |r| {r[0], r[1]} }
    intv_ovlp_idx(ivs, idx, st, en).map { |i| a[i] }
  end

  def self.cmd_sveval(args : Array(String)) : Int32
    min_flt = 30; min_size = 50; max_size = 100_000; win_size = 500
    print_err = false; print_match = false; bed_fn : String? = nil
    len_diff_ratio = 0.5
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-f"          ; i += 1; min_flt = parse_num(args[i])
      when "-i"          ; i += 1; min_size = parse_num(args[i])
      when "-x"          ; i += 1; max_size = parse_num(args[i])
      when "-w"          ; i += 1; win_size = parse_num(args[i])
      when "-d"          ; i += 1; len_diff_ratio = args[i].to_f
      when "-r"          ; i += 1; bed_fn = args[i]
      when "-e"          ; print_err = true
      when "-p"          ; print_match = true
      when "-h", "--help"; STDERR.puts "Usage: paftools sveval [options] <base.vcf> <call.vcf>"; return 0
      else                 rest << args[i]
      end
      i += 1
    end
    if rest.size < 2
      STDERR.puts "Usage: paftools sveval [options] <base.vcf> <call.vcf>"
      STDERR.puts "Options:"
      STDERR.puts "  -r FILE    confident region in BED []"
      STDERR.puts "  -f INT     min length to discard [#{min_flt}]"
      STDERR.puts "  -i INT     min SV length [#{min_size}]"
      STDERR.puts "  -x INT     max SV length [#{max_size}]"
      STDERR.puts "  -w INT     fuzzy windown size [#{win_size}]"
      STDERR.puts "  -d FLOAT   max allele diff if there is a single allele in the window [#{len_diff_ratio}]"
      STDERR.puts "  -e         print errors"
      return 1
    end

    bed : Hash(String, Array({Int32, Int32}))? = nil
    bed_idx : Hash(String, Array(Int32))? = nil
    if fn = bed_fn
      bed, bed_idx = sveval_read_bed(fn)
    end

    compare_vcf = ->(v0 : Hash(String, Array(SvRec)), v1 : Hash(String, Array(SvRec)), label : String) {
      m = 0; n = 0
      v1.each do |x, a1|
        a0 = v0[x]?
        a0_idx = a0 ? sveval_idx(a0) : nil
        a1.each do |rec|
          next if rec[3] < min_size
          n += 1
          next unless a0 && a0_idx
          ws = win_size + (rec[3] >> 1)
          st = rec[0] > ws ? rec[0] - ws : 0
          b = sveval_ovlp(a0, a0_idx, st, rec[1] + ws)
          n_ins = 0; n_del = 0; sv_del : Int32? = nil; sv_ins : Int32? = nil
          b.each do |brec|
            bsvlen = brec[2]
            if bsvlen && bsvlen < 0
              n_del += 1; sv_del = -bsvlen
            elsif bsvlen && bsvlen > 0
              n_ins += 1; sv_ins = bsvlen
            end
            puts ["MA", x, rec[0], rec[1], rec[2], brec[0], brec[1], brec[2]].join('\t') if print_match
          end
          match = false
          rsvlen = rec[2]
          if rsvlen && rsvlen > 0 # insertion
            if n_ins == 1
              diff = (sv_ins.not_nil! - rec[3]).abs
              match = true if diff < min_size || diff.to_f / rec[3] < len_diff_ratio
            elsif n_ins > 1
              match = true # multiple insertions; ambiguous
            end
          elsif rsvlen && rsvlen < 0
            if n_del == 1
              diff = (sv_del.not_nil! - rec[3]).abs
              match = true if diff < min_size || diff.to_f / rec[3] < len_diff_ratio
            elsif n_del > 1
              match = true # multiple deletions; ambiguous
            end
          end
          if match
            m += 1
          elsif print_err
            if (rsvlen && rsvlen > 0 && n_ins > 0) || (rsvlen && rsvlen < 0 && n_del > 0)
              puts ["MM", x, rec[0], rec[1], rec[2]].join('\t')
            end
            puts [label, x, rec[0], rec[1], rec[2]].join('\t')
          end
        end
      end
      {n, m}
    }

    v_base = sveval_read_vcf(rest[0], bed, bed_idx, min_flt, max_size)
    v_call = sveval_read_vcf(rest[1], bed, bed_idx, min_flt, max_size)
    fn_stat = compare_vcf.call(v_call, v_base, "FN")
    fp_stat = compare_vcf.call(v_base, v_call, "FP")
    # NOTE: k8's print() joins multiple arguments with TAB, not space.
    # Uses js_fixed (not raw "%.6f" %) because these ratios can be 0/0 when
    # no SV in the callset/baseline reaches min_size: JS's toFixed prints
    # "NaN", but Crystal's sprintf prints the different string "nan".
    puts ["SN", fn_stat[0], fn_stat[1], js_fixed(fn_stat[1].to_f / fn_stat[0], 6)].join('\t')
    puts ["PC", fp_stat[0], fp_stat[1], js_fixed(fp_stat[1].to_f / fp_stat[0], 6)].join('\t')
    puts ["F1", js_fixed((fn_stat[1].to_f / fn_stat[0] + fp_stat[1].to_f / fp_stat[0]) / 2, 6)].join('\t')
    0
  end
end
