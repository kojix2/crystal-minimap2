module Minimap2
  # ---------------------------------------------------------------------------
  # Seed collection and filtering — port of seed.c
  # ---------------------------------------------------------------------------

  # Filter overly-frequent minimizers from a query minimizer set.
  # Mirrors mm_seed_mz_flt().
  def self.seed_mz_flt(mv : Array(Mm128), q_occ_max : Int32, q_occ_frac : Float32) : Nil
    return if mv.size <= q_occ_max || q_occ_frac <= 0.0_f32 || q_occ_max <= 0

    # Sort a copy by minimizer hash, track original indices
    a = mv.map_with_index { |m, i| Mm128.new(m.x, i.to_u64) }
    radix_sort_128x(a)

    # Mark overly frequent minimizers (zero out their x)
    st = 0
    i = 1
    while i <= a.size
      if i == a.size || a[i].x != a[st].x
        cnt = i - st
        if cnt > q_occ_max && cnt > mv.size * q_occ_frac
          (st...i).each { |j| mv[a[j].y] = Mm128.new(0_u64, mv[a[j].y].y) }
        end
        st = i
      end
      i += 1
    end

    # Compact (remove zeroed entries)
    j = 0
    mv.each_with_index do |m, k|
      mv[j] = mv[k]
      j += 1 if m.x != 0
    end
    mv.delete_at(j, mv.size - (j)) if mv.size > (j)
  end

  # Collect all seed matches between query minimizers and the index.
  # Mirrors mm_seed_collect_all().
  # Returns an Array(MmSeed) with all matches (including high-freq ones).
  def self.seed_collect_all(mi : MmIdx, mv : Array(Mm128)) : Array(MmSeed)
    seeds = [] of MmSeed
    mv.each_with_index do |p, i|
      q_pos = p.y.to_u32
      q_span = (p.x & 0xff).to_u32
      cr = mi.get(p.x >> 8)
      next unless cr
      next if cr.empty?
      seed = MmSeed.new(
        n: cr.size,
        q_pos: q_pos,
        q_span: q_span,
        cr: cr.dup,
        seg_id: (p.y >> 32).to_u32
      )
      # mark tandem if adjacent minimizers share the same hash
      if i > 0 && p.x >> 8 == mv[i - 1].x >> 8
        seed = MmSeed.new(seed.n, seed.q_pos, seed.q_span, seed.cr, seed.seg_id, false, true)
      end
      if i < mv.size - 1 && p.x >> 8 == mv[i + 1].x >> 8
        seed = MmSeed.new(seed.n, seed.q_pos, seed.q_span, seed.cr, seed.seg_id, false, true)
      end
      seeds << seed
    end
    seeds
  end

  MAX_MAX_HIGH_OCC = 128

  # Select up to max_high_occ high-frequency minimizers per streak.
  # Mirrors mm_seed_select().
  def self.seed_select(seeds : Array(MmSeed), qlen : Int32,
                       max_occ : Int32, max_max_occ : Int32, dist : Int32) : Nil
    n = seeds.size
    return if n <= 1

    n_hi = seeds.count { |s| s.n > max_occ }
    return if n_hi == 0

    last0 = -1
    i = 0
    while i <= n
      if i == n || seeds[i].n <= max_occ
        if i - last0 > 1
          ps = last0 < 0 ? 0 : (seeds[last0].q_pos >> 1).to_i32
          pe = i == n ? qlen : (seeds[i].q_pos >> 1).to_i32
          st = last0 + 1
          en = i
          max_high_occ = ((pe - ps).to_f / dist + 0.499).to_i32
          if max_high_occ > 0
            max_high_occ = MAX_MAX_HIGH_OCC if max_high_occ > MAX_MAX_HIGH_OCC
            # Pick top max_high_occ by n (using a heap or just sort)
            run = (st...en).to_a
            # Sort run by n descending; take first max_high_occ
            run.sort_by! { |j| -seeds[j].n }
            chosen = run.first(max_high_occ).to_set
            (st...en).each do |j|
              # invert: flag those NOT chosen
              if chosen.includes?(j)
                seeds[j] = MmSeed.new(seeds[j].n, seeds[j].q_pos, seeds[j].q_span, seeds[j].cr, seeds[j].seg_id, false, seeds[j].is_tandem?)
              else
                seeds[j] = MmSeed.new(seeds[j].n, seeds[j].q_pos, seeds[j].q_span, seeds[j].cr, seeds[j].seg_id, true, seeds[j].is_tandem?)
              end
            end
          end
          # Filter anything exceeding max_max_occ
          (st...en).each do |j|
            if seeds[j].n > max_max_occ
              seeds[j] = MmSeed.new(seeds[j].n, seeds[j].q_pos, seeds[j].q_span, seeds[j].cr, seeds[j].seg_id, true, seeds[j].is_tandem?)
            end
          end
        end
        last0 = i
      end
      i += 1
    end
  end

  # Full seed collection pipeline (mirrors mm_collect_matches).
  # Returns:
  #   seeds    — filtered MmSeed array
  #   n_a      — total reference anchor count
  #   rep_len  — repetitive region length in query
  #   mini_pos — (q_span<<32 | q_pos) for each non-filtered seed
  def self.collect_matches(
    qlen : Int32, max_occ : Int32, max_max_occ : Int32, dist : Int32,
    mi : MmIdx, mv : Array(Mm128),
  ) : {Array(MmSeed), Int64, Int32, Array(UInt64)}
    seeds = seed_collect_all(mi, mv)

    if dist > 0 && max_max_occ > max_occ
      seed_select(seeds, qlen, max_occ, max_max_occ, dist)
    else
      seeds.each_with_index do |s, i|
        if s.n > max_occ
          seeds[i] = MmSeed.new(s.n, s.q_pos, s.q_span, s.cr, s.seg_id, true, s.is_tandem?)
        end
      end
    end

    n_a = 0_i64
    rep_len = 0_i32
    rep_st = 0_i32
    rep_en = 0_i32
    mini_pos = [] of UInt64
    out = [] of MmSeed

    seeds.each do |q|
      if (@@dbg_flag & DBG_SEED_FREQ) != 0
        STDERR.printf("SF\t%d\t%d\t%d\n", q.q_pos >> 1, q.n, q.flt? ? 1 : 0)
      end
      if q.flt?
        en = (q.q_pos >> 1) + 1
        st = en - q.q_span.to_i32
        if st > rep_en
          rep_len += rep_en - rep_st
          rep_st = st; rep_en = en
        else
          rep_en = en
        end
      else
        n_a += q.n
        mini_pos << (q.q_span.to_u64 << 32 | (q.q_pos >> 1).to_u64)
        out << q
      end
    end
    rep_len += rep_en - rep_st

    {out, n_a, rep_len, mini_pos}
  end
end
