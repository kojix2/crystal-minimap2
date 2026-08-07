module Paftools
  # ── Numeric formatting ──────────────────────────────────────────────────
  # Mirrors JS Number#toFixed(2), including its NaN/Infinity string output
  # (Crystal's "%.2f" has no such special-casing).
  def self.js_fixed2(v : Float64) : String
    js_fixed(v, 2)
  end

  def self.js_fixed(v : Float64, digits : Int32) : String
    return "NaN" if v.nan?
    return (v > 0 ? "Infinity" : "-Infinity") if v.infinite?
    "%.#{digits}f" % v
  end

  # ── I/O ──────────────────────────────────────────────────────────────────

  def self.open_in(fn : String, &)
    if fn == "-"
      yield STDIN
    elsif fn.ends_with?(".gz") || fn.ends_with?(".bgz")
      File.open(fn, "rb") { |file| Compress::Gzip::Reader.open(file) { |gzip| yield gzip } }
    else
      File.open(fn) { |file| yield file }
    end
  end

  # ── Interval operations ───────────────────────────────────────────────────
  # Sorted {start,end} pairs with a link-back index for O(n) overlap queries
  # (mirrors Interval.sort / merge / index_end / find_ovlp in paftools.js).

  def self.intv_sort(a : Array({Int32, Int32}))
    a.sort_by! { |intv| {intv[0], intv[1]} }
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

  # ── Numeric argument parsing (mirrors JS paf_parseNum) ────────────────────
  # Parses sizes like "1m", "500k", "2.5g" used by -l/-g/-f/-i/-x/-w options.

  def self.parse_num(s : String) : Int32
    m = /^(\d*\.?\d*)([mMgGkK]?)/.match(s)
    return 0 unless m
    x = m[1].to_f
    case m[2]
    when "k", "K" then x *= 1_000
    when "m", "M" then x *= 1_000_000
    when "g", "G" then x *= 1_000_000_000
    end
    (x + 0.499).floor.to_i32
  end

  # ── VCF ALT-allele length helper (mirrors JS _paf_get_alen) ───────────────
  # Returns {alen, min_abs_diff, max_abs_diff} for a VCF record's ALT alleles,
  # preferring INFO/SVLEN when present, else falling back to the length delta
  # between REF and each non-symbolic ALT allele.
  def self.paf_get_alen(t : Array(String)) : {Int32?, Int32, Int32}
    svlen : Int32? = nil
    if m = /(^|;)SVLEN=(-?\d+)/.match(t[7])
      svlen = m[2].to_i
    end
    s = t[4].split(',')
    min_abs_diff = 1 << 30
    max_abs_diff = 0
    alen : Int32? = nil
    if (sv = svlen) && sv != 0
      alen = sv
      min_abs_diff = max_abs_diff = sv > 0 ? sv : -sv
    end
    rlen = t[3].size
    s.each do |allele|
      next if /^<\S+>$/.matches?(allele)
      diff = allele.size - rlen
      abs_diff = diff > 0 ? diff : -diff
      min_abs_diff = [min_abs_diff, abs_diff].min
      if max_abs_diff < abs_diff
        max_abs_diff = abs_diff
        alen = diff
      end
    end
    {alen, min_abs_diff, max_abs_diff}
  end

  # ── Merge overlapping regions and return covered length ───────────────────

  def self.cov_len(regs : Array({Int32, Int32})) : Int64
    return 0_i64 if regs.empty?
    s = regs.sort_by { |reg| reg[0] }
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

  def self.read_fasta(fn : String) : {Hash(String, String), Array({String, Int32})}
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
end
