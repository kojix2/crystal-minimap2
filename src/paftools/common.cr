module Paftools
  # ── I/O ──────────────────────────────────────────────────────────────────

  def self.open_in(fn : String, &)
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

  def self.cov_len(regs : Array({Int32, Int32})) : Int64
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
