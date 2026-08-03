module Minimap2
  # ---------------------------------------------------------------------------
  # Global state  (mirrors misc.c globals)
  # ---------------------------------------------------------------------------
  @@verbose : Int32 = 1 # 0=silent, 1=error, 2=warning, 3=message
  @@dbg_flag : Int32 = 0
  @@t0 : Time::Instant = Time.instant # program-start monotonic timestamp

  def self.verbose : Int32
    @@verbose
  end

  def self.verbose=(v : Int32)
    @@verbose = v
  end

  def self.dbg_flag : Int32
    @@dbg_flag
  end

  def self.dbg_flag=(v : Int32)
    @@dbg_flag = v
  end

  # Seconds elapsed since program start (used as realtime0 baseline: always 0.0).
  def self.realtime0 : Float64
    0.0
  end

  # ---------------------------------------------------------------------------
  # Truncate UInt64 to lower 32 bits and reinterpret as Int32.
  # Equivalent to the C cast (int32_t)v — truncates without raising.
  # ---------------------------------------------------------------------------
  @[AlwaysInline]
  def self.u64_to_i32(v : UInt64) : Int32
    (v & 0xffffffff_u64).to_u32.unsafe_as(Int32)
  end

  # ---------------------------------------------------------------------------
  # Seconds elapsed since program start on the monotonic clock
  # (equivalent to realtime() in C — used only for logging, not wall time).
  # ---------------------------------------------------------------------------
  def self.realtime : Float64
    Time.instant.duration_since(@@t0).total_seconds
  end

  # ---------------------------------------------------------------------------
  # CPU time consumed by the current process (user + system), in seconds
  # ---------------------------------------------------------------------------
  def self.cputime : Float64
    t = Process.times
    t.utime + t.stime
  end

  # ---------------------------------------------------------------------------
  # Peak resident-set size in bytes (Linux only; 0 on other platforms)
  # ---------------------------------------------------------------------------
  def self.peakrss : Int64
    {% if flag?(:linux) %}
      ru = LibC::Rusage.new
      LibC.getrusage(LibC::RUSAGE_SELF, pointerof(ru))
      ru.ru_maxrss.to_i64 * 1024
    {% else %}
      0_i64
    {% end %}
  end

  # ---------------------------------------------------------------------------
  # Utility: round x up to the next power of two (32-bit version of
  # kroundup32 from ksort.h / mmpriv.h)
  # ---------------------------------------------------------------------------
  def self.kroundup32(x : Int32) : Int32
    x -= 1
    x |= x >> 1
    x |= x >> 2
    x |= x >> 4
    x |= x >> 8
    x |= x >> 16
    x + 1
  end

  def self.kroundup64(x : Int64) : Int64
    x -= 1
    x |= x >> 1
    x |= x >> 2
    x |= x >> 4
    x |= x >> 8
    x |= x >> 16
    x |= x >> 32
    x + 1
  end

  # ---------------------------------------------------------------------------
  # Fast log₂ approximation (mg_log2 from mmpriv.h)
  # Only works for x >= 2.
  # ---------------------------------------------------------------------------
  def self.mg_log2(x : Float32) : Float32
    # bit-trick: extract biased exponent and fractional part
    bits = x.unsafe_as(UInt32)
    log_2 = ((bits >> 23) & 0xff_u32).to_f32 - 128.0_f32
    bits &= ~(0xff_u32 << 23)
    bits |= 127_u32 << 23
    z = bits.unsafe_as(Float32)
    log_2 + (-0.34484843_f32 * z + 2.02466578_f32) * z - 0.67487759_f32
  end

  # ---------------------------------------------------------------------------
  # Radix sort for Array(Mm128) by key = x (8-bit radix, 8 passes).
  # Uses a pre-allocated class-level buffer to avoid per-call GC allocation.
  # ---------------------------------------------------------------------------
  @@sort_buf = [] of Mm128

  def self.radix_sort_128x(a : Array(Mm128)) : Nil
    radix_sort_128x_impl(a, 0, a.size)
  end

  def self.radix_sort_128x(a : Array(Mm128), from : Int32, to : Int32) : Nil
    radix_sort_128x_impl(a, from, to)
  end

  private def self.radix_sort_128x_impl(a : Array(Mm128), from : Int32, to : Int32) : Nil
    n = to - from
    return if n <= 1

    # Ensure class-level buffer is large enough (never shrinks).
    if @@sort_buf.size < n
      @@sort_buf = Array(Mm128).new(n, Mm128.max)
    end
    buf = @@sort_buf

    count = StaticArray(Int32, 256).new(0)

    8.times do |pass|
      shift = pass * 8
      count.fill(0)

      # Count frequencies
      from.upto(to - 1) do |i|
        count[((a[i].x >> shift) & 0xff).to_i32] += 1
      end

      # Skip if all in one bucket
      non_zero = 0
      256.times { |bi| non_zero += 1 if count[bi] > 0 }
      next if non_zero <= 1

      # Prefix sum
      sum = 0
      256.times do |bi|
        c = count[bi]
        count[bi] = sum
        sum += c
      end

      # Scatter to buffer
      from.upto(to - 1) do |i|
        b = ((a[i].x >> shift) & 0xff).to_i32
        buf[count[b]] = a[i]
        count[b] += 1
      end

      # Copy back
      from.upto(to - 1) do |i|
        a[i] = buf[i - from]
      end
    end
  end

  # ---------------------------------------------------------------------------
  # Radix sort for Array(UInt64) — 8-bit radix, 8 passes.
  # ---------------------------------------------------------------------------
  @@sort_buf64 = [] of UInt64

  def self.radix_sort_64(a : Array(UInt64)) : Nil
    radix_sort_64_impl(a, 0, a.size)
  end

  def self.radix_sort_64(a : Array(UInt64), from : Int32, to : Int32) : Nil
    radix_sort_64_impl(a, from, to)
  end

  private def self.radix_sort_64_impl(a : Array(UInt64), from : Int32, to : Int32) : Nil
    n = to - from
    return if n <= 1

    if @@sort_buf64.size < n
      @@sort_buf64 = Array(UInt64).new(n, 0_u64)
    end
    buf = @@sort_buf64

    count = StaticArray(Int32, 256).new(0)

    8.times do |pass|
      shift = pass * 8
      count.fill(0)

      from.upto(to - 1) do |i|
        count[((a[i] >> shift) & 0xff).to_i32] += 1
      end

      non_zero = 0
      256.times { |bi| non_zero += 1 if count[bi] > 0 }
      next if non_zero <= 1

      sum = 0
      256.times do |bi|
        c = count[bi]
        count[bi] = sum
        sum += c
      end

      from.upto(to - 1) do |i|
        b = ((a[i] >> shift) & 0xff).to_i32
        buf[count[b]] = a[i]
        count[b] += 1
      end

      from.upto(to - 1) do |i|
        a[i] = buf[i - from]
      end
    end
  end

  # ---------------------------------------------------------------------------
  # k-th smallest element (in-place sort, no .dup needed for callers)
  # ---------------------------------------------------------------------------
  def self.ks_ksmall_uint32(a : Array(UInt32), kk : Int32) : UInt32
    a.sort!
    a[kk]
  end
end
