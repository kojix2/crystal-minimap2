module Minimap2
  # ---------------------------------------------------------------------------
  # Symmetric DUST — low-complexity region masking
  # Port of sdust.c
  #
  # Returns an array of [start, end) intervals of masked regions,
  # packed as UInt64 (start<<32 | end).
  # ---------------------------------------------------------------------------

  private SD_WLEN = 3
  private SD_WTOT = 1 << (SD_WLEN * 2)  # 64
  private SD_WMSK = SD_WTOT - 1          # 63

  # A "perfect interval" found during the DUST scan.
  private record PerfIntv, start : Int32, finish : Int32, r : Int32, l : Int32

  # ---------------------------------------------------------------------------
  # Reusable per-call buffer
  # ---------------------------------------------------------------------------
  class SdustBuf
    property w : Array(Int32)    # sliding window queue (acts as deque)
    property big_p : Array(PerfIntv)  # list of perfect intervals
    property res : Array(UInt64)  # masked regions result

    def initialize
      @w   = [] of Int32
      @big_p = [] of PerfIntv
      @res = [] of UInt64
    end

    def reset
      @w.clear
      @big_p.clear
      @res.clear
    end
  end

  # ---------------------------------------------------------------------------
  # shift_window: update the sliding window and scoring variables.
  # Mirrors the static shift_window() in sdust.c.
  # ---------------------------------------------------------------------------
  private def self.sdust_shift_window(
    t : Int32, w : Array(Int32), big_t : Int32, big_w : Int32,
    l : Int32, rw : Int32, rv : Int32, cw : Array(Int32), cv : Array(Int32)
  ) : {Int32, Int32, Int32}
    if w.size >= big_w - SD_WLEN + 1
      s = w.shift  # pop from front
      rw -= (cw[s] -= 1)
      if l > w.size + 1   # +1 because we already shifted
        l -= 1
        rv -= (cv[s] -= 1)
      end
    end
    w << t
    l += 1
    rw += (cw[t] += 1) - 1  # rw += old cw[t] (before increment)
    rv += (cv[t] += 1) - 1
    # ensure cv[t]*10 <= T*2 (T*2 because we're using integer arithmetic)
    while cv[t] * 10 > big_t * 2
      s = w[w.size - l]
      rv -= (cv[s] -= 1)
      l -= 1
      break if s == t  # stop once we've found t in the L-window
    end
    {l, rw, rv}
  end

  # ---------------------------------------------------------------------------
  # save_masked_regions: append any perfect intervals that have fallen out
  # of the current window to the result list, merging overlapping intervals.
  # ---------------------------------------------------------------------------
  private def self.sdust_save_masked(res : Array(UInt64), big_p : Array(PerfIntv), start : Int32) : Nil
    return if big_p.empty? || big_p.last.start >= start
    p = big_p.last
    if res.size > 0
      s = (res.last >> 32).to_i32
      f = res.last.to_u32.to_i32
      if p.start <= f
        res[-1] = (s.to_u64 << 32) | (f > p.finish ? f : p.finish).to_u64
        # remove trailing perfect intervals that are now out of window
        while big_p.size > 0 && big_p.last.start < start
          big_p.pop
        end
        return
      end
    end
    res << (p.start.to_u64 << 32 | p.finish.to_u64)
    # remove trailing intervals that fell out of window
    while big_p.size > 0 && big_p.last.start < start
      big_p.pop
    end
  end

  # ---------------------------------------------------------------------------
  # find_perfect: look for new perfect intervals in the window.
  # Mirrors find_perfect() in sdust.c.
  # ---------------------------------------------------------------------------
  private def self.sdust_find_perfect(
    big_p : Array(PerfIntv), w : Array(Int32),
    big_t : Int32, start : Int32, big_l : Int32, rv : Int32, cv : Array(Int32)
  ) : Nil
    c  = cv.dup  # local copy of the cv counts
    r  = rv
    max_r = 0
    max_l = 0
    # Scan from oldest to newest in L-window (in reverse order in the array)
    i = w.size - big_l - 1
    while i >= 0
      t     = w[i]
      r    += c[t]
      c[t] += 1
      new_r = r
      new_l = w.size - i - 1
      if new_r * 10 > big_t * new_l
        # Find insertion point (big_p is sorted by descending start)
        j = 0
        while j < big_p.size && big_p[j].start >= i + start
          pi = big_p[j]
          max_r, max_l = pi.r, pi.l if max_r == 0 || pi.r * max_l > max_r * pi.l
          j += 1
        end
        if max_r == 0 || new_r * max_l >= max_r * new_l
          max_r = new_r; max_l = new_l
          big_p.insert(j, PerfIntv.new(
            start:  i + start,
            finish: w.size + (SD_WLEN - 1) + start,
            r:      new_r,
            l:      new_l
          ))
        end
      end
      i -= 1
    end
  end

  # ---------------------------------------------------------------------------
  # sdust_core: run DUST on one sequence, using a reusable buffer.
  # Returns a Slice view into buf.res.
  # ---------------------------------------------------------------------------
  def self.sdust_core(seq : Bytes, l_seq : Int32, big_t : Int32, big_w : Int32,
                      buf : SdustBuf) : Array(UInt64)
    buf.reset
    cw = Array(Int32).new(SD_WTOT, 0)
    cv = Array(Int32).new(SD_WTOT, 0)
    rv  = 0
    rw  = 0
    big_l = 0
    t     = 0
    l     = 0  # length of current run of unambiguous bases

    i = 0
    while i <= l_seq
      b = i < l_seq ? SEQ_NT4_TABLE[seq[i].to_i].to_i : 4
      if b < 4
        l += 1
        t = ((t << 2) | b) & SD_WMSK
        if l >= SD_WLEN
          start = ([l - big_w, 0].max) + (i + 1 - l)
          sdust_save_masked(buf.res, buf.big_p, start)
          big_l, rw, rv = sdust_shift_window(t, buf.w, big_t, big_w, big_l, rw, rv, cw, cv)
          if rw * 10 > big_l * big_t
            sdust_find_perfect(buf.big_p, buf.w, big_t, start, big_l, rv, cv)
          end
        end
      else
        start = ([l - big_w + 1, 0].max) + (i + 1 - l)
        while buf.big_p.size > 0
          sdust_save_masked(buf.res, buf.big_p, start)
          start += 1
        end
        l = t = 0
        cw.fill(0); cv.fill(0)
        big_l = rw = rv = 0
        buf.w.clear
      end
      i += 1
    end

    buf.res
  end

  # ---------------------------------------------------------------------------
  # sdust: allocate a buffer, run sdust_core, return results.
  # Equivalent to sdust() in sdust.c.
  # ---------------------------------------------------------------------------
  def self.sdust(seq : Bytes, l_seq : Int32 = -1, big_t : Int32 = 20, big_w : Int32 = 64) : Array(UInt64)
    len = l_seq < 0 ? seq.size : l_seq
    buf = SdustBuf.new
    sdust_core(seq, len, big_t, big_w, buf).dup
  end
end
