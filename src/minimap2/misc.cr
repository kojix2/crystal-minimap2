module Minimap2
  # ---------------------------------------------------------------------------
  # Global state  (mirrors misc.c globals)
  # ---------------------------------------------------------------------------
  @@verbose    : Int32         = 1   # 0=silent, 1=error, 2=warning, 3=message
  @@dbg_flag   : Int32         = 0
  @@t0         : Time::Instant = Time.instant  # program-start monotonic timestamp

  def self.verbose : Int32;       @@verbose;         end
  def self.verbose=(v : Int32);   @@verbose = v;     end
  def self.dbg_flag : Int32;      @@dbg_flag;        end
  def self.dbg_flag=(v : Int32);  @@dbg_flag = v;    end

  # Seconds elapsed since program start (used as realtime0 baseline: always 0.0).
  def self.realtime0 : Float64;   0.0;               end

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
  # Radix sort for Array(Mm128) by key = x (equivalent to radix_sort_128x)
  # In Crystal we just use Array#sort_by! which is introsort; the radix sort
  # matters for performance but not correctness.
  # ---------------------------------------------------------------------------
  def self.radix_sort_128x(a : Array(Mm128)) : Nil
    a.sort_by!(&.x)
  end

  def self.radix_sort_128x(a : Array(Mm128), from : Int32, to : Int32) : Nil
    # sort sub-slice [from, to)
    slice = a[from...to]
    slice.sort_by!(&.x)
    a[from...to] = slice
  end

  # ---------------------------------------------------------------------------
  # Radix sort for Array(UInt64) (equivalent to radix_sort_64)
  # ---------------------------------------------------------------------------
  def self.radix_sort_64(a : Array(UInt64)) : Nil
    a.sort!
  end

  def self.radix_sort_64(a : Array(UInt64), from : Int32, to : Int32) : Nil
    slice = a[from...to]
    slice.sort!
    a[from...to] = slice
  end

  # ---------------------------------------------------------------------------
  # k-th smallest element (equivalent to ks_ksmall_uint32_t)
  # Uses a partial sort / selection algorithm.
  # ---------------------------------------------------------------------------
  def self.ks_ksmall_uint32(a : Array(UInt32), kk : Int32) : UInt32
    arr = a.dup
    arr.sort!
    arr[kk]
  end
end
