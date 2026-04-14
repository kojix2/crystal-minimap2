module Minimap2
  # ---------------------------------------------------------------------------
  # Linear chaining — port of lchain.c
  # ---------------------------------------------------------------------------

  INT32_MIN = -2147483648_i32
  INT32_MAX =  2147483647_i32

  # ---------------------------------------------------------------------------
  # comput_sc: per-anchor scoring function (mirrors static comput_sc())
  # ---------------------------------------------------------------------------
  private def self.comput_sc(
    ai : Mm128, aj : Mm128,
    max_dist_x : Int32, max_dist_y : Int32, bw : Int32,
    chn_pen_gap : Float32, chn_pen_skip : Float32,
    is_cdna : Bool, n_seg : Int32,
  ) : Int32
    dq = u64_to_i32(ai.y) - u64_to_i32(aj.y)
    return INT32_MIN if dq <= 0 || dq > max_dist_x

    dr = u64_to_i32(ai.x - aj.x)
    sidi = ((ai.y & SEED_SEG_MASK) >> SEED_SEG_SHIFT).to_i32
    sidj = ((aj.y & SEED_SEG_MASK) >> SEED_SEG_SHIFT).to_i32

    return INT32_MIN if sidi == sidj && (dr == 0 || dq > max_dist_y)

    dd = dr > dq ? dr - dq : dq - dr
    return INT32_MIN if sidi == sidj && dd > bw
    return INT32_MIN if n_seg > 1 && !is_cdna && sidi == sidj && dr > max_dist_y

    dg = dr < dq ? dr : dq
    q_span = ((aj.y >> 32) & 0xff).to_i32
    sc = q_span < dg ? q_span : dg

    if dd != 0 || dg > q_span
      lin_pen = chn_pen_gap * dd.to_f32 + chn_pen_skip * dg.to_f32
      log_pen = dd >= 1 ? mg_log2((dd + 1).to_f32) : 0.0_f32

      if is_cdna || sidi != sidj
        if sidi != sidj && dr == 0
          sc += 1
        elsif dr > dq || sidi != sidj
          sc -= [lin_pen, log_pen].min.to_i32
        else
          sc -= (lin_pen + 0.5_f32 * log_pen).to_i32
        end
      else
        sc -= (lin_pen + 0.5_f32 * log_pen).to_i32
      end
    end
    sc
  end

  # ---------------------------------------------------------------------------
  # comput_sc_simple: simplified scoring (for RMQ chaining, mirrors comput_sc_simple)
  # ---------------------------------------------------------------------------
  private def self.comput_sc_simple(
    ai : Mm128, aj : Mm128,
    chn_pen_gap : Float32, chn_pen_skip : Float32,
    exact : Pointer(Bool)?, width : Pointer(Int32),
  ) : Int32
    dq = u64_to_i32(ai.y) - u64_to_i32(aj.y)
    dr = u64_to_i32(ai.x - aj.x)
    dd = dr > dq ? dr - dq : dq - dr
    dg = dr < dq ? dr : dq
    q_span = ((aj.y >> 32) & 0xff).to_i32
    sc = q_span < dg ? q_span : dg
    width.value = dd
    exact.try(&.value=(((dd == 0) && dg <= q_span)))

    if dd != 0 || dq > q_span
      lin_pen = chn_pen_gap * dd.to_f32 + chn_pen_skip * dg.to_f32
      log_pen = dd >= 1 ? mg_log2((dd + 1).to_f32) : 0.0_f32
      sc -= (lin_pen + 0.5_f32 * log_pen).to_i32
    end
    sc
  end

  # ---------------------------------------------------------------------------
  # mg_chain_bk_end: walk back from anchor k to find end of chain region.
  # Mirrors the static mg_chain_bk_end() in lchain.c.
  # ---------------------------------------------------------------------------
  private def self.chain_bk_end(
    max_drop : Int32, z : Array(Mm128),
    f : Array(Int32), p : Array(Int64),
    t : Array(Int32), k : Int64,
  ) : Int64
    i = z[k].y.to_i64
    end_i = -1_i64
    max_i = i
    max_s = 0_i32

    return i if i < 0 || t[i] != 0

    loop do
      t[i] = 2
      end_i = i
      i = p[i]
      s = i < 0 ? u64_to_i32(z[k].x) : (u64_to_i32(z[k].x) - f[i])
      if s > max_s
        max_s = s
        max_i = i
      elsif max_s - s > max_drop
        break
      end
      break if i < 0 || t[i] != 0
    end

    # Reset modified t[]
    i = z[k].y.to_i64
    while i >= 0 && i != end_i
      t[i] = 0
      i = p[i]
    end
    max_i
  end

  # ---------------------------------------------------------------------------
  # mg_chain_backtrack: extract chains from the DP arrays.
  # Returns (u, v) or (nil, nil) if no chains found.
  # u[i] = score<<32 | count;  v[] = anchor indices (reversed)
  # ---------------------------------------------------------------------------
  def self.chain_backtrack(
    n : Int64, f : Array(Int32), p : Array(Int64),
    v : Array(Int32), t : Array(Int32),
    min_cnt : Int32, min_sc : Int32, max_drop : Int32,
  ) : {Array(UInt64), Int32, Array(Int32), Int32}
    n_z = 0_i64
    0.upto(n - 1) { |i| n_z += 1 if f[i] >= min_sc }

    return {[] of UInt64, 0, [] of Int32, 0} if n_z == 0

    # Build z[] = sorted-by-score (score, original_index)
    z = [] of Mm128
    0.upto(n - 1) do |i|
      z << Mm128.new(f[i].unsafe_as(UInt32).to_u64, i.to_u64) if f[i] >= min_sc
    end
    radix_sort_128x(z)

    t.fill(0)
    n_v = 0
    n_u = 0

    # First pass: count
    (n_z - 1).downto(0) do |k|
      next unless t[z[k].y] == 0
      n_v0 = n_v
      end_i = chain_bk_end(max_drop, z, f, p, t, k)
      i = z[k].y.to_i64
      while i != end_i
        n_v += 1
        t[i] = 1
        i = p[i]
      end
      sc = i < 0 ? u64_to_i32(z[k].x) : u64_to_i32(z[k].x) - f[i]
      if sc >= min_sc && n_v > n_v0 && n_v - n_v0 >= min_cnt
        n_u += 1
      else
        n_v = n_v0
      end
    end

    u_arr = Array(UInt64).new(n_u, 0_u64)
    v_arr = Array(Int32).new(n_v, 0)
    t.fill(0)
    n_v = 0
    n_u = 0

    # Second pass: populate
    (n_z - 1).downto(0) do |k|
      next unless t[z[k].y] == 0
      n_v0 = n_v
      end_i = chain_bk_end(max_drop, z, f, p, t, k)
      i = z[k].y.to_i64
      while i != end_i
        v_arr[n_v] = i.to_i32
        n_v += 1
        t[i] = 1
        i = p[i]
      end
      sc = i < 0 ? u64_to_i32(z[k].x) : u64_to_i32(z[k].x) - f[i]
      if sc >= min_sc && n_v > n_v0 && n_v - n_v0 >= min_cnt
        u_arr[n_u] = sc.to_u64 << 32 | (n_v - n_v0).to_u64
        n_u += 1
      else
        n_v = n_v0
      end
    end

    {u_arr, n_u, v_arr, n_v}
  end

  # Rearrange anchors from v[] into contiguous runs ordered by chain, then
  # sort chains by target position (mirrors compact_a).
  private def self.compact_a(
    n_u : Int32, u : Array(UInt64),
    n_v : Int32, v : Array(Int32), a : Array(Mm128),
  ) : Array(Mm128)
    b = Array(Mm128).new(n_v, Mm128.max)
    k = 0
    0.upto(n_u - 1) do |i|
      k0 = k
      ni = u64_to_i32(u[i]) # lower 32 bits = count
      0.upto(ni - 1) { |j| b[k] = a[v[k0 + (ni - j - 1)]]; k += 1 }
    end

    # Sort chains by target position of first anchor
    w = Array(Mm128).new(n_u) { Mm128.max }
    pos = 0
    0.upto(n_u - 1) do |i|
      w[i] = Mm128.new(b[pos].x, (pos.to_u64 << 32) | i.to_u64)
      pos += u64_to_i32(u[i])
    end
    radix_sort_128x(w)

    u2 = Array(UInt64).new(n_u, 0_u64)
    result = Array(Mm128).new(n_v, Mm128.max)
    k = 0
    0.upto(n_u - 1) do |i|
      j = u64_to_i32(w[i].y) # original chain index (lower 32 bits)
      ni = u64_to_i32(u[j])
      src = (w[i].y >> 32).to_i32 # position in b[]
      u2[i] = u[j]
      ni.times do |off|
        result[k] = b[src + off]
        k += 1
      end
    end
    u[0...n_u] = u2

    result
  end

  # ---------------------------------------------------------------------------
  # mg_lchain_dp: DP-based linear chaining (mirrors mg_lchain_dp in lchain.c).
  #
  # Input:
  #   a[].x = rev<<63 | tid<<32 | tpos
  #   a[].y = flags<<40 | q_span<<32 | q_pos
  # Output:
  #   n_u (via return value): number of chains
  #   u[i] = score<<32 | #anchors
  #   returns rearranged a[] (may be empty if no chains)
  # ---------------------------------------------------------------------------
  def self.lchain_dp(
    max_dist_x : Int32, max_dist_y : Int32, bw : Int32,
    max_skip : Int32, max_iter : Int32,
    min_cnt : Int32, min_sc : Int32,
    chn_pen_gap : Float32, chn_pen_skip : Float32,
    is_cdna : Bool, n_seg : Int32,
    a : Array(Mm128), u_out : Array(UInt64),
  ) : Array(Mm128)
    n = a.size.to_i64
    u_out.clear

    return [] of Mm128 if n == 0

    max_dist_x = bw if max_dist_x < bw
    max_dist_y = bw if max_dist_y < bw && !is_cdna
    max_drop = is_cdna ? INT32_MAX : bw

    f = Array(Int32).new(n, 0)
    p = Array(Int64).new(n, -1_i64)
    v = Array(Int32).new(n, 0)
    t = Array(Int32).new(n, 0)

    st = 0_i64
    max_ii = -1_i64

    n.times do |i|
      ai = a[i]
      max_j = -1_i64
      max_f = ((ai.y >> 32) & 0xff).to_i32
      n_skip = 0

      # Advance st to discard out-of-range anchors
      while st < i && (a[i].x >> 32 != a[st].x >> 32 || a[i].x > a[st].x + max_dist_x)
        st += 1
      end
      st = i - max_iter if i - st > max_iter

      # DP inner loop
      end_j = st
      j = i - 1
      while j >= st
        sc = comput_sc(ai, a[j], max_dist_x, max_dist_y, bw,
          chn_pen_gap, chn_pen_skip, is_cdna, n_seg)
        if sc != INT32_MIN
          sc += f[j]
          if sc > max_f
            max_f = sc; max_j = j
            n_skip -= 1 if n_skip > 0
          elsif t[j] == i.to_i32
            n_skip += 1
            break if n_skip > max_skip
          end
          t[p[j]] = i.to_i32 if p[j] >= 0
        end
        end_j = j
        j -= 1
      end

      # Check max_ii anchor
      if max_ii < 0 || a[i].x - a[max_ii].x > max_dist_x
        max = INT32_MIN
        max_ii = -1_i64
        (st...i).each do |anc_j|
          if max < f[anc_j]
            max = f[anc_j]; max_ii = anc_j
          end
        end
      end
      if max_ii >= 0 && max_ii < end_j
        tmp = comput_sc(ai, a[max_ii], max_dist_x, max_dist_y, bw,
          chn_pen_gap, chn_pen_skip, is_cdna, n_seg)
        if tmp != INT32_MIN && max_f < tmp + f[max_ii]
          max_f = tmp + f[max_ii]; max_j = max_ii
        end
      end

      f[i] = max_f; p[i] = max_j
      v[i] = max_j >= 0 && v[max_j] > max_f ? v[max_j] : max_f
      max_ii = i if max_ii < 0 || (a[i].x - a[max_ii].x <= max_dist_x && f[max_ii] < f[i])
    end

    u_arr, n_u, v_arr, n_v = chain_backtrack(n, f, p, v, t, min_cnt, min_sc, max_drop)
    u_arr[0...n_u].each { |u_val| u_out << u_val }

    return [] of Mm128 if n_u == 0

    compact_a(n_u, u_out, n_v, v_arr, a)
  end

  # ---------------------------------------------------------------------------
  # mg_lchain_rmq: RMQ-based linear chaining (mirrors mg_lchain_rmq in lchain.c).
  #
  # For Crystal we implement the same algorithm using a sorted structure
  # instead of KRMQ (which is a C macro-based RMQ tree).
  # We use a simple sorted array with binary search for the range-max query.
  # ---------------------------------------------------------------------------

  # A comparable element for the RMQ tree
  private record LcElem, y : Int32, i : Int64, pri : Float64

  # Naive RMQ over a sorted list: find element with max priority (min -pri)
  # in the range where y is within [lo_y, hi_y].
  private def self.rmq_query(tree : Array(LcElem), lo_y : Int32, hi_y : Int32) : LcElem?
    best : LcElem? = nil
    tree.each do |e|
      next if e.y < lo_y || e.y > hi_y
      if best.nil? || e.pri < best.pri # lower pri = higher f score
        best = e
      end
    end
    best
  end

  def self.lchain_rmq(
    max_dist : Int32, max_dist_inner : Int32, bw : Int32,
    max_chn_skip : Int32, cap_rmq_size : Int32,
    min_cnt : Int32, min_sc : Int32,
    chn_pen_gap : Float32, chn_pen_skip : Float32,
    a : Array(Mm128), u_out : Array(UInt64),
  ) : Array(Mm128)
    n = a.size.to_i64
    u_out.clear

    return [] of Mm128 if n == 0

    max_dist = bw if max_dist < bw
    max_dist_inner = 0 if max_dist_inner < 0
    max_dist_inner = max_dist if max_dist_inner > max_dist
    max_drop = bw

    f = Array(Int32).new(n, 0)
    p = Array(Int64).new(n, -1_i64)
    t = Array(Int32).new(n, 0)
    v = Array(Int32).new(n, 0)

    root = [] of LcElem # sorted by y
    root_inner = [] of LcElem

    st = 0_i64
    st_inner = 0_i64
    i0 = 0_i64

    n.times do |i|
      ai = a[i]
      max_j = -1_i64
      q_span = ((ai.y >> 32) & 0xff).to_i32
      max_f = q_span

      # Add in-range anchors to RMQ tree
      if i0 < i && a[i0].x != ai.x
        (i0...i).each do |j|
          elem = LcElem.new(
            y: u64_to_i32(a[j].y),
            i: j,
            pri: -(f[j].to_f64 + 0.5 * chn_pen_gap * (u64_to_i32(a[j].x) + u64_to_i32(a[j].y)))
          )
          # Insert sorted by y
          ins = root.bsearch_index { |e| e.y >= elem.y } || root.size
          root.insert(ins, elem)
          if max_dist_inner > 0
            ins2 = root_inner.bsearch_index { |e| e.y >= elem.y } || root_inner.size
            root_inner.insert(ins2, elem)
          end
        end
        i0 = i
      end

      # Evict out-of-range anchors from root
      while st < i && (ai.x >> 32 != a[st].x >> 32 ||
            ai.x > a[st].x + max_dist ||
            root.size > cap_rmq_size)
        root.reject! { |e| e.i == st }
        st += 1
      end
      if max_dist_inner > 0
        while st_inner < i && (ai.x >> 32 != a[st_inner].x >> 32 ||
              ai.x > a[st_inner].x + max_dist_inner ||
              root_inner.size > cap_rmq_size)
          root_inner.reject! { |e| e.i == st_inner }
          st_inner += 1
        end
      end

      # RMQ query on main tree
      lo_y = u64_to_i32(ai.y) - max_dist
      hi_y = u64_to_i32(ai.y)

      if q = rmq_query(root, lo_y, hi_y)
        j = q.i
        w = 0
        sc = f[j] + comput_sc_simple(ai, a[j], chn_pen_gap, chn_pen_skip, nil, pointerof(w))
        if w <= bw && sc > max_f
          max_f = sc; max_j = j
        end

        # Inner tree query for exact matches
        if max_dist_inner > 0 && u64_to_i32(ai.y) > 0
          lo2 = u64_to_i32(ai.y) - max_dist_inner
          n_skip = 0

          root_inner.reverse_each do |elem|
            next if elem.y >= u64_to_i32(ai.y) # must be strictly less
            break if elem.y < lo2

            jj = elem.i
            sc2 = f[jj] + comput_sc_simple(ai, a[jj], chn_pen_gap, chn_pen_skip, nil, pointerof(w))
            if w <= bw
              if sc2 > max_f
                max_f = sc2; max_j = jj
                n_skip -= 1 if n_skip > 0
              elsif t[jj] == i.to_i32
                n_skip += 1
                break if n_skip > max_chn_skip
              end
              t[p[jj]] = i.to_i32 if p[jj] >= 0
            end
          end
        end
      end

      f[i] = max_f; p[i] = max_j
      v[i] = max_j >= 0 && v[max_j] > max_f ? v[max_j] : max_f
    end

    u_arr, n_u, v_arr, n_v = chain_backtrack(n, f, p, v, t, min_cnt, min_sc, max_drop)
    u_arr[0...n_u].each { |u_val| u_out << u_val }

    return [] of Mm128 if n_u == 0

    compact_a(n_u, u_out, n_v, v_arr, a)
  end
end
