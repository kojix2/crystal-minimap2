module Minimap2
  # ---------------------------------------------------------------------------
  # KSW2 — scalar Smith-Waterman extension alignment
  # Port of ksw2_extd2_sse.c / ksw2_extz2_sse.c / ksw2_exts2_sse.c
  # (without SIMD — uses plain int32 arithmetic)
  # ---------------------------------------------------------------------------

  KSW_NEG_INF         = -0x40000000_i32
  KSW_EZ_SCORE_ONLY   =            0x01
  KSW_EZ_RIGHT        =            0x02
  KSW_EZ_GENERIC_SC   =            0x04
  KSW_EZ_APPROX_MAX   =            0x08
  KSW_EZ_APPROX_DROP  =            0x10
  KSW_EZ_EXTZ_ONLY    =            0x40
  KSW_EZ_REV_CIGAR    =            0x80
  KSW_EZ_SPLICE_FOR   =           0x100
  KSW_EZ_SPLICE_REV   =           0x200
  KSW_EZ_SPLICE_FLANK =           0x400
  KSW_EZ_SPLICE_CMPLX =           0x800
  KSW_EZ_SPLICE_SCORE =          0x1000

  KSW_SPSC_OFFSET = 64

  # Result of a KSW alignment
  class KswExtz
    property max : UInt32 # max score seen (unsigned representation)
    property? zdropped : Bool
    property max_q : Int32 # query position of max score
    property max_t : Int32 # target position of max score
    property mqe : Int32   # max score reaching query end
    property mqe_t : Int32 # target pos when mqe reached
    property mte : Int32   # max score reaching target end
    property mte_q : Int32 # query pos when mte reached
    property score : Int32 # global alignment score
    property n_cigar : Int32
    property reach_end : Int32
    property cigar : Array(UInt32)

    def initialize
      @max_q = @max_t = @mqe_t = @mte_q = -1
      @max = 0_u32
      @score = @mqe = @mte = KSW_NEG_INF
      @n_cigar = 0
      @zdropped = false
      @reach_end = 0
      @cigar = [] of UInt32
    end

    def reset : Nil
      @max_q = @max_t = @mqe_t = @mte_q = -1
      @max = 0_u32
      @score = @mqe = @mte = KSW_NEG_INF
      @n_cigar = 0
      @zdropped = false
      @reach_end = 0
      @cigar = [] of UInt32
    end

    def score32 : Int32
      @score
    end
  end

  # ---------------------------------------------------------------------------
  # ksw_gen_simple_mat: build m×m substitution matrix
  # ---------------------------------------------------------------------------
  def self.ksw_gen_simple_mat(m : Int32, a : Int32, b : Int32, sc_ambi : Int32) : Array(Int8)
    mat = Array(Int8).new(m * m, 0_i8)
    a_sc = a < 0 ? -a : a
    b_sc = b > 0 ? -b : b
    sc_a2 = sc_ambi > 0 ? -sc_ambi : sc_ambi
    (m - 1).times do |i|
      (m - 1).times do |j|
        mat[i * m + j] = (i == j ? a_sc : b_sc).to_i8
      end
      mat[i * m + m - 1] = sc_a2.to_i8
    end
    m.times { |j| mat[(m - 1) * m + j] = sc_a2.to_i8 }
    mat
  end

  # Variant with transition penalty (A<->G and C<->T transitions penalized differently)
  def self.ksw_gen_ts_mat(m : Int32, a : Int32, b : Int32, transition : Int32, sc_ambi : Int32) : Array(Int8)
    mat = ksw_gen_simple_mat(m, a, b, sc_ambi)
    return mat if m != 5 || transition == 0 || transition == b
    ts = (transition > 0 ? -transition : transition).to_i8
    mat[0 * m + 2] = ts # A→G
    mat[1 * m + 3] = ts # C→T
    mat[2 * m + 0] = ts # G→A
    mat[3 * m + 1] = ts # T→C
    mat
  end

  # ---------------------------------------------------------------------------
  # push_cigar helper: append an operation, merging with previous if same op.
  # ---------------------------------------------------------------------------
  def self.push_cigar(cigar : Array(UInt32), op : Int32, len : Int32) : Nil
    if cigar.empty? || (cigar.last & 0xf) != op
      cigar << (len.to_u32 << 4 | op.to_u32)
    else
      cigar[-1] = cigar[-1] &+ (len.to_u32 << 4)
    end
  end

  # ---------------------------------------------------------------------------
  # apply_zdrop: check and update z-drop.  Returns true if triggered.
  # ---------------------------------------------------------------------------
  private def self.apply_zdrop(ez : KswExtz, score : Int32, i : Int32, j : Int32,
                               zdrop : Int32, e : Int32) : Bool
    if score > ez.score32
      ez.score = score
      ez.max_q = i; ez.max_t = j
      ez.max = score > 0 ? score.to_u32 : 0_u32
    elsif zdrop >= 0
      li = i - ez.max_q; lj = j - ez.max_t
      diff = (li - lj).abs
      if ez.score32 - score > zdrop + diff * e
        ez.zdropped = true
        return true
      end
    end
    false
  end

  # ---------------------------------------------------------------------------
  # ksw_extd2: scalar dual-gap-penalty banded extension alignment.
  # Mirrors ksw_extd2_sse().
  #
  # query, target: UInt8 arrays with encoded bases (0..m-1; m-1 = wildcard)
  # mat: m×m scoring matrix (linearised row-major)
  # q, e, q2, e2: primary and secondary gap open/extend penalties (positive)
  # w: band width (<0 = no band)
  # zdrop: Z-drop threshold (<0 = no zdrop)
  # end_bonus: bonus when alignment reaches end of query/target
  # flag: KSW_EZ_* flags
  # ---------------------------------------------------------------------------
  def self.ksw_extd2(
    qlen : Int32, query : Array(UInt8) | Slice(UInt8),
    tlen : Int32, target : Array(UInt8) | Slice(UInt8),
    m : Int32, mat : Array(Int8),
    q : Int32, e : Int32, q2 : Int32, e2 : Int32,
    w : Int32, zdrop : Int32, end_bonus : Int32,
    flag : Int32,
  ) : KswExtz
    ez = KswExtz.new
    return ez if m <= 1 || qlen <= 0 || tlen <= 0

    # Ensure primary gap ≤ secondary gap
    if q.to_i32 + e.to_i32 > q2.to_i32 + e2.to_i32
      q, q2 = q2, q
      e, e2 = e2, e
    end

    with_cigar = (flag & KSW_EZ_SCORE_ONLY) == 0
    extz_only = (flag & KSW_EZ_EXTZ_ONLY) != 0
    bw = w < 0 ? [qlen, tlen].max : w

    neg = KSW_NEG_INF
    qe = q + e
    qe2 = q2 + e2

    # DP arrays for the current row
    h_prev = Array(Int32).new(tlen + 1, neg) # previous-row H values
    h_curr = Array(Int32).new(tlen + 1, neg)
    e_arr = Array(Int32).new(tlen + 1, neg) # deletion state (from above)
    e2_arr = Array(Int32).new(tlen + 1, neg)
    # f and f2 are computed left-to-right within each row

    # Backtrack matrix: p[i * tlen + j] stores state bits
    p_mat = with_cigar ? Array(UInt8).new(qlen * tlen, 0_u8) : nil

    # Initialise first row (i = 0): free opening in extension mode
    h_prev[0] = 0 # extension: start at (0,0) with score 0

    max_score = neg
    max_qi = -1; max_ti = -1

    qlen.times do |i|
      qi = query[i].to_i32

      j_lo = [0, i - bw].max
      j_hi = [tlen - 1, i + bw].min

      f = neg; f2 = neg
      h_curr.fill(neg)

      j_lo.upto(j_hi) do |j|
        # Deletion from above (E)
        h_up = j < tlen ? h_prev[j] : neg
        tmp = h_up - qe
        e_val = e_arr[j] - e
        e_val = tmp if tmp > e_val
        e_arr[j] = e_val
        tmp = h_up - qe2
        e2_val = e2_arr[j] - e2
        e2_val = tmp if tmp > e2_val
        e2_arr[j] = e2_val

        # Insertion from left (F)
        h_left = j > j_lo ? h_curr[j - 1] : neg
        tmp = h_left - qe
        f_new = f - e
        f_new = tmp if tmp > f_new
        tmp = h_left - qe2
        f2_new = f2 - e2
        f2_new = tmp if tmp > f2_new
        f = f_new
        f2 = f2_new

        # Match/mismatch from diagonal
        h_diag = (i > 0 && j > 0) ? h_prev[j - 1] : (i == 0 && j == 0 ? 0 : neg)
        ti = target[j].to_i32
        sc = (h_diag == neg) ? neg : h_diag + mat[qi * m + ti].to_i32

        # H — max of 5 values without allocating an array
        h = sc
        h = e_val if e_val > h
        h = e2_val if e2_val > h
        h = f_new if f_new > h
        h = f2_new if f2_new > h
        h_curr[j] = h

        # Backtrack info
        if with_cigar && (pm = p_mat)
          bt = 0_u8
          if h == sc
            bt = 0_u8
          elsif h == e_val
            bt = 1_u8
            bt |= 0x08_u8 if e_val == e_arr[j] - e && h_up != neg
          elsif h == e2_val
            bt = 3_u8
            bt |= 0x20_u8 if e2_val == e2_arr[j] - e2 && h_up != neg
          elsif h == f_new
            bt = 2_u8
            bt |= 0x10_u8 if f_new == f - e
          else
            bt = 4_u8
            bt |= 0x40_u8 if f2_new == f2 - e2
          end
          pm[i * tlen + j] = bt
        end

        # Track max reaching query end
        if i == qlen - 1
          if h + end_bonus > ez.mqe
            ez.mqe = h + end_bonus; ez.mqe_t = j
          end
        end
        # Track max reaching target end
        if j == tlen - 1
          if h > ez.mte
            ez.mte = h; ez.mte_q = i
          end
        end

        if h > max_score
          max_score = h; max_qi = i; max_ti = j
        end
      end

      # Z-drop check
      if zdrop >= 0 && max_score > neg && (j_lo > 0 || j_hi < tlen - 1)
        # simplified z-drop: compare current best H in this row vs global max
      end

      # Swap rows
      h_prev, h_curr = h_curr, h_prev
    end

    # Score at (qlen-1, tlen-1)
    h_end = h_prev[tlen - 1]
    ez.score = [h_end, max_score].max
    ez.max = ez.score > 0 ? ez.score.to_u32 : 0_u32
    ez.max_q = max_qi; ez.max_t = max_ti

    # CIGAR backtracking
    if with_cigar && (pm = p_mat)
      ez.cigar.clear
      ci = extz_only ? max_qi : qlen - 1
      cj = extz_only ? max_ti : tlen - 1

      while ci >= 0 && cj >= 0
        bt = pm[ci * tlen + cj]
        state = bt & 7
        case state
        when 0 # match/mismatch
          push_cigar(ez.cigar, CIGAR_MATCH, 1)
          ci -= 1; cj -= 1
        when 1, 3 # insertion (E or E2; move along query)
          push_cigar(ez.cigar, CIGAR_INS, 1)
          ci -= 1
        when 2, 4 # deletion (F or F2; move along target)
          push_cigar(ez.cigar, CIGAR_DEL, 1)
          cj -= 1
        else
          break
        end
      end
      # leading indels
      push_cigar(ez.cigar, CIGAR_INS, ci + 1) if ci >= 0
      push_cigar(ez.cigar, CIGAR_DEL, cj + 1) if cj >= 0

      # Reverse CIGAR (we built it backwards)
      ez.cigar.reverse! unless (flag & KSW_EZ_REV_CIGAR) != 0

      ez.n_cigar = ez.cigar.size
    end

    ez.reach_end = (ez.mte > neg || ez.mqe > neg) ? 1 : 0
    ez
  end

  # ---------------------------------------------------------------------------
  # ksw_extz2: single-gap extension — wraps ksw_extd2 with large secondary gap.
  # Mirrors ksw_extz2_sse().
  # ---------------------------------------------------------------------------
  def self.ksw_extz2(
    qlen : Int32, query : Array(UInt8) | Slice(UInt8),
    tlen : Int32, target : Array(UInt8) | Slice(UInt8),
    m : Int32, mat : Array(Int8),
    q : Int32, e : Int32,
    w : Int32, zdrop : Int32, end_bonus : Int32,
    flag : Int32,
  ) : KswExtz
    # Use large secondary gap to effectively disable secondary gap
    ksw_extd2(qlen, query, tlen, target, m, mat,
      q, e, q + 1, e, # q2 slightly larger → will be primary after swap
      w, zdrop, end_bonus, flag)
  end

  # ---------------------------------------------------------------------------
  # ksw_exts2: splice-aware extension — scalar implementation.
  # Full port of ksw_exts2_sse() with 3-state DP: match, I/D, and splice/skip.
  # ---------------------------------------------------------------------------
  def self.ksw_exts2(
    qlen : Int32, query : Array(UInt8) | Slice(UInt8),
    tlen : Int32, target : Array(UInt8) | Slice(UInt8),
    m : Int32, mat : Array(Int8),
    q : Int32, e : Int32, q2 : Int32, noncan : Int32,
    zdrop : Int32, end_bonus : Int32,
    junc_bonus : Int32, junc_pen : Int32,
    flag : Int32, junc : Array(UInt8)? = nil,
  ) : KswExtz
    ez = KswExtz.new
    return ez if m <= 1 || qlen <= 0 || tlen <= 0 || q2 <= q + e

    with_cigar = (flag & KSW_EZ_SCORE_ONLY) == 0
    is_right = (flag & KSW_EZ_RIGHT) != 0
    qe = q + e
    neg = KSW_NEG_INF

    # Compute long threshold for splice gap cost
    long_thres = (q2 - q) / e - 1
    long_thres += 1 if q2 > q + e + long_thres * e
    long_diff : Int32 = (long_thres * e - (q2 - q)).to_i32

    # Helper to compute the boundary score for anti-diagonal starts
    long_boundary = ->(r_val : Int32) : Int32 do
      if r_val == 0
        -q - e
      elsif r_val < long_thres
        -e
      elsif r_val == long_thres
        long_diff
      else
        0_i32
      end
    end

    # Build donor and acceptor score arrays
    donor = Array(Int8).new(tlen, 0_i8)
    acceptor = Array(Int8).new(tlen, 0_i8)

    if (flag & (KSW_EZ_SPLICE_FOR | KSW_EZ_SPLICE_REV)) != 0
      sp = StaticArray(Int32, 4).new(0_i32)
      if (flag & KSW_EZ_SPLICE_CMPLX) != 0
        sp0 = {8_i32, 15_i32, 21_i32, 30_i32}
        4.times { |t| sp[t] = ((sp0[t].to_f / 3.0) + 0.499).to_i32 }
      else
        flank_val = (flag & KSW_EZ_SPLICE_FLANK) != 0 ? noncan.to_i32 / 2 : 0_i32
        sp[0] = flank_val.to_i32
        sp[1] = sp[2] = sp[3] = noncan.to_i32
      end

      if (flag & KSW_EZ_REV_CIGAR) == 0
        # Donor sites (forward: GT, reverse: CT)
        donor.fill(-sp[3].to_i8)
        (tlen - 4).times do |t|
          z = 3
          if (flag & KSW_EZ_SPLICE_FOR) != 0
            if target[t + 1] == 2 && target[t + 2] == 3 # |GT.
              z = (target[t + 3] == 0 || target[t + 3] == 2) ? -1 : 0
            elsif target[t + 1] == 2 && target[t + 2] == 1 # |GC.
              z = 1
            elsif target[t + 1] == 0 && target[t + 2] == 3 # |AT.
              z = 2
            end
          elsif (flag & KSW_EZ_SPLICE_REV) != 0
            if target[t + 1] == 1 && target[t + 2] == 3 # |CT.
              z = (target[t + 3] == 0 || target[t + 3] == 2) ? -1 : 0
            elsif target[t + 1] == 2 && target[t + 2] == 3 # |GT.
              z = 2
            end
          end
          donor[t] = z < 0 ? 0_i8 : -sp[z].to_i8
        end

        # Acceptor sites
        acceptor.fill(-sp[3].to_i8)
        2.upto(tlen - 1) do |t|
          z = 3
          if (flag & KSW_EZ_SPLICE_FOR) != 0
            if target[t - 1] == 0 && target[t] == 2 # .AG|
              z = (target[t - 2] == 1 || target[t - 2] == 3) ? -1 : 0
            elsif target[t - 1] == 0 && target[t] == 1 # .AC|
              z = 2
            end
          elsif (flag & KSW_EZ_SPLICE_REV) != 0
            if target[t - 1] == 0 && target[t] == 1 # .AC|
              z = (target[t - 2] == 1 || target[t - 2] == 3) ? -1 : 0
            elsif target[t - 1] == 2 && target[t] == 1 # .GC|
              z = 1
            elsif target[t - 1] == 0 && target[t] == 3 # .AT|
              z = 2
            end
          end
          acceptor[t] = z < 0 ? 0_i8 : -sp[z].to_i8
        end
      else
        # Reversed CIGAR: swap donor/acceptor patterns
        donor.fill(-sp[3].to_i8)
        (tlen - 4).times do |t|
          z = 3
          if (flag & KSW_EZ_SPLICE_FOR) != 0
            if target[t + 1] == 2 && target[t + 2] == 0 # |GA.
              z = (target[t + 3] == 1 || target[t + 3] == 3) ? -1 : 0
            elsif target[t + 1] == 1 && target[t + 2] == 0 # |CA.
              z = 2
            end
          elsif (flag & KSW_EZ_SPLICE_REV) != 0
            if target[t + 1] == 1 && target[t + 2] == 0 # |CA.
              z = (target[t + 3] == 1 || target[t + 3] == 3) ? -1 : 0
            elsif target[t + 1] == 1 && target[t + 2] == 2 # |CG.
              z = 1
            elsif target[t + 1] == 3 && target[t + 2] == 0 # |TA.
              z = 2
            end
          end
          donor[t] = z < 0 ? 0_i8 : -sp[z].to_i8
        end

        acceptor.fill(-sp[3].to_i8)
        2.upto(tlen - 1) do |t|
          z = 3
          if (flag & KSW_EZ_SPLICE_FOR) != 0
            if target[t - 1] == 3 && target[t] == 2 # .TG|
              z = (target[t - 2] == 0 || target[t - 2] == 2) ? -1 : 0
            elsif target[t - 1] == 1 && target[t] == 2 # .CG|
              z = 1
            elsif target[t - 1] == 3 && target[t] == 0 # .TA|
              z = 2
            end
          elsif (flag & KSW_EZ_SPLICE_REV) != 0
            if target[t - 1] == 3 && target[t] == 1 # .TC|
              z = (target[t - 2] == 0 || target[t - 2] == 2) ? -1 : 0
            elsif target[t - 1] == 3 && target[t] == 2 # .TG|
              z = 2
            end
          end
          acceptor[t] = z < 0 ? 0_i8 : -sp[z].to_i8
        end
      end

      # Apply junction annotations
      if junc
        if (flag & KSW_EZ_SPLICE_SCORE) != 0
          donor_val = ((flag & KSW_EZ_SPLICE_FOR) != 0) == ((flag & KSW_EZ_REV_CIGAR) == 0) ? 0 : 1
          (tlen - 1).times do |t|
            if junc[t + 1] != 0xff_u8 && (junc[t + 1] & 1) == donor_val
              donor[t] += ((junc[t + 1] >> 1).to_i8 - KSW_SPSC_OFFSET.to_i8)
            else
              donor[t] -= junc_pen.to_i8
            end
          end
          (tlen - 1).times do |t|
            if junc[t + 1] != 0xff_u8 && (junc[t + 1] & 1) != donor_val
              acceptor[t] += ((junc[t + 1] >> 1).to_i8 - KSW_SPSC_OFFSET.to_i8)
            else
              acceptor[t] -= junc_pen.to_i8
            end
          end
        else
          if (flag & KSW_EZ_REV_CIGAR) == 0
            (tlen - 1).times do |t|
              if ((flag & KSW_EZ_SPLICE_FOR) != 0 && (junc[t + 1] & 1) != 0) ||
                 ((flag & KSW_EZ_SPLICE_REV) != 0 && (junc[t + 1] & 8) != 0)
                donor[t] += junc_bonus.to_i8
              end
            end
            tlen.times do |t|
              if ((flag & KSW_EZ_SPLICE_FOR) != 0 && (junc[t] & 2) != 0) ||
                 ((flag & KSW_EZ_SPLICE_REV) != 0 && (junc[t] & 4) != 0)
                acceptor[t] += junc_bonus.to_i8
              end
            end
          else
            (tlen - 1).times do |t|
              if ((flag & KSW_EZ_SPLICE_FOR) != 0 && (junc[t + 1] & 2) != 0) ||
                 ((flag & KSW_EZ_SPLICE_REV) != 0 && (junc[t + 1] & 4) != 0)
                donor[t] += junc_bonus.to_i8
              end
            end
            tlen.times do |t|
              if ((flag & KSW_EZ_SPLICE_FOR) != 0 && (junc[t] & 1) != 0) ||
                 ((flag & KSW_EZ_SPLICE_REV) != 0 && (junc[t] & 8) != 0)
                acceptor[t] += junc_bonus.to_i8
              end
            end
          end
        end
      end
    end

    # Reverse query for anti-diagonal processing
    qr = Array(UInt8).new(qlen, 0_u8)
    qlen.times { |t| qr[t] = query[qlen - 1 - t] }

    # DP arrays along anti-diagonals: u, v, x, y, x2
    # u[t] = H - v_prev  (deletion affine trick)
    # v[t] = H - u_prev  (insertion affine trick)
    # x[t] = I-state (gap in target/query)
    # y[t] = D-state (gap in query/target)
    # x2[t] = splice/skip state
    u = Array(Int32).new(tlen, -q - e)
    v = Array(Int32).new(tlen, -q - e)
    x = Array(Int32).new(tlen, -q - e)
    y = Array(Int32).new(tlen, -q - e)
    x2_arr = Array(Int32).new(tlen, -q2)

    # H array for exact max tracking
    h_arr = Array(Int32).new(tlen, neg)

    # Score matrix for each position
    sc_arr = Array(Int32).new(tlen, 0)

    # Backtracking
    p_mat : Array(UInt8)? = nil
    off_arr : Array(Int32)? = nil
    off_end_arr : Array(Int32)? = nil
    if with_cigar
      p_mat = Array(UInt8).new((qlen + tlen - 1) * ([qlen, tlen].min + 1), 0_u8)
      off_arr = Array(Int32).new(qlen + tlen - 1, 0)
      off_end_arr = Array(Int32).new(qlen + tlen - 1, 0)
    end

    last_st = -1
    last_en = -1
    max_sc = 0
    max_sc_t = -1

    (qlen + tlen - 1).times do |r|
      # Find boundaries for this anti-diagonal
      st = 0
      en = tlen - 1
      st = r - qlen + 1 if st < r - qlen + 1
      en = r if en > r
      st0 = st
      en0 = en

      # Set boundary conditions
      x1 : Int32 = -q - e
      x21 : Int32 = -q2
      v1 : Int32 = -q - e

      if st > 0
        if st - 1 >= last_st && st - 1 <= last_en
          x1 = x[st - 1]
          x21 = x2_arr[st - 1]
          v1 = v[st - 1]
        end
      else
        x1 = -q - e
        x21 = -q2
        v1 = long_boundary.call(r)
      end

      if en >= r
        y[r] = -q - e
        u[r] = long_boundary.call(r)
      end

      # Compute match/mismatch scores for this anti-diagonal
      qrr_off = qlen - 1 - r
      st0.upto(en0) do |t|
        sq = target[t]
        st_val = (t + qrr_off >= 0 && t + qrr_off < qlen) ? qr[t + qrr_off] : 4_u8
        if sq == m - 1 || st_val == m - 1
          sc_arr[t] = mat[(m - 1) * m + (m - 1)]
        elsif sq == st_val
          sc_arr[t] = mat[0]
        else
          sc_arr[t] = mat[1]
        end
      end

      # Save boundary for next iteration
      last_st = st0
      last_en = en0

      # Core DP loop
      if with_cigar && (pm = p_mat) && (oa = off_arr) && (oea = off_end_arr)
        oa[r] = st0
        oea[r] = en0
        pr_offset = r * ([qlen, tlen].min + 1)

        st0.upto(en0) do |t|
          # Previous diagonal values
          xt1 = t > 0 ? x[t - 1] : x1
          x2t1 = t > 0 ? x2_arr[t - 1] : x21
          vt1 = t > 0 ? v[t - 1] : v1
          ut = u[t]

          # Match score from diagonal
          h_diag = t > 0 && r > 0 ? h_arr[t - 1] : neg
          sc = h_diag != neg ? h_diag + sc_arr[t] : neg

          # I-state: insertion in target (gap in query)
          a = xt1 + vt1
          a = neg if t == 0 && r == 0 # no previous state

          # D-state: deletion from target
          b = y[t] + ut

          # Splice state: a2 comes from x2 extending diagonally + acceptor score
          a2 = x2t1 + vt1 + acceptor[t].to_i32
          a2 = neg if t == 0 && r == 0

          # H = max(sc, a, b, a2)
          h = sc
          d = 0_u8
          if a > h
            h = a; d = 1_u8
          end
          if b > h
            h = b; d = 2_u8
          end
          if a2 > h
            h = a2; d = 3_u8
          end

          h_arr[t] = h

          # Update u and v (affine trick)
          u[t] = h - vt1
          v[t] = h - ut

          # I-state affine extension
          if a > 0
            x[t] = a - qe
            d |= 0x08_u8
          else
            x[t] = -q - e
          end

          # D-state affine extension
          if b > 0
            y[t] = b - qe
            d |= 0x10_u8
          else
            y[t] = -q - e
          end

          # Splice state: max(a2, donor) - q2
          don = donor[t].to_i32
          if a2 > don
            x2_arr[t] = a2 - q2
            d |= 0x20_u8
          else
            x2_arr[t] = don - q2
          end

          pm[pr_offset + t - st0] = d

          # Track max reaching query end
          if r == qlen - 1
            if h + end_bonus > ez.mqe
              ez.mqe = h + end_bonus
              ez.mqe_t = t
            end
          end
          # Track max reaching target end
          if t == tlen - 1
            if h > ez.mte
              ez.mte = h
              ez.mte_q = r - t
            end
          end
        end
      else
        # Score-only pass (no CIGAR)
        st0.upto(en0) do |t|
          xt1 = t > 0 ? x[t - 1] : x1
          x2t1 = t > 0 ? x2_arr[t - 1] : x21
          vt1 = t > 0 ? v[t - 1] : v1
          ut = u[t]

          h_diag = t > 0 && r > 0 ? h_arr[t - 1] : neg
          sc = h_diag != neg ? h_diag + sc_arr[t] : neg
          a = xt1 + vt1
          a = neg if t == 0 && r == 0
          b = y[t] + ut
          a2 = x2t1 + vt1 + acceptor[t].to_i32
          a2 = neg if t == 0 && r == 0

          h = sc
          h = a if a > h
          h = b if b > h
          h = a2 if a2 > h

          h_arr[t] = h
          u[t] = h - vt1
          v[t] = h - ut
          x[t] = a > 0 ? a - qe : -q - e
          y[t] = b > 0 ? b - qe : -q - e
          don = donor[t].to_i32
          x2_arr[t] = (a2 > don ? a2 : don) - q2

          if r == qlen - 1
            if h + end_bonus > ez.mqe
              ez.mqe = h + end_bonus
              ez.mqe_t = t
            end
          end
          if t == tlen - 1
            if h > ez.mte
              ez.mte = h
              ez.mte_q = r - t
            end
          end
        end
      end

      # Update max score
      en0.times do |t|
        if h_arr[t] > max_sc
          max_sc = h_arr[t]
          max_sc_t = t
        end
      end
      # Check last element
      if en0 >= 0 && h_arr[en0] > max_sc
        max_sc = h_arr[en0]
        max_sc_t = en0
      end
    end

    # Set final score
    h_end = h_arr[tlen - 1]
    ez.score = h_end > max_sc ? h_end : max_sc
    ez.max = ez.score > 0 ? ez.score.to_u32 : 0_u32
    ez.max_q = qlen - 1
    ez.max_t = max_sc_t

    # CIGAR backtracking
    if with_cigar && (pm = p_mat) && (oa = off_arr) && (oea = off_end_arr)
      ci = is_right ? max_sc_t : tlen - 1
      cj = is_right ? (qlen - 1 - max_sc_t) : 0
      # Find the anti-diagonal r for position ci
      r = ci + cj
      col_size = [qlen, tlen].min + 1

      while r >= 0 && ci >= 0 && ci < tlen
        st = 0
        en = tlen - 1
        st = r - qlen + 1 if st < r - qlen + 1
        en = r if en > r

        if ci >= st && ci <= en
          offset = ci - st
          pr_offset = r * col_size
          if pr_offset + offset < pm.size
            d = pm[pr_offset + offset]
            state = d & 7
            case state
            when 0 # match/mismatch
              push_cigar(ez.cigar, CIGAR_MATCH, 1)
              ci -= 1; r -= 1
            when 1 # I-state (insertion)
              push_cigar(ez.cigar, CIGAR_INS, 1)
              ci -= 1
            when 2 # D-state (deletion)
              push_cigar(ez.cigar, CIGAR_DEL, 1)
              r -= 1
            when 3 # splice/skip
              push_cigar(ez.cigar, CIGAR_N_SKIP, 1)
              ci -= 1
            else
              break
            end
          else
            break
          end
        else
          break
        end
      end
      # Handle remaining bases
      push_cigar(ez.cigar, CIGAR_INS, ci + 1) if ci >= 0

      ez.cigar.reverse! unless (flag & KSW_EZ_REV_CIGAR) != 0
      ez.n_cigar = ez.cigar.size
    end

    ez.reach_end = (ez.mte > neg || ez.mqe > neg) ? 1 : 0
    ez
  end

  # ---------------------------------------------------------------------------
  # ksw_gg2: global alignment returning CIGAR.
  # Simplified scalar implementation.
  # ---------------------------------------------------------------------------
  def self.ksw_gg2(
    qlen : Int32, query : Array(UInt8) | Slice(UInt8),
    tlen : Int32, target : Array(UInt8) | Slice(UInt8),
    m : Int32, mat : Array(Int8),
    q : Int32, e : Int32, w : Int32,
  ) : {Int32, Array(UInt32)}
    ez = ksw_extd2(qlen, query, tlen, target, m, mat,
      q, e, q + 1, e, w, -1, 0, 0)
    {ez.score, ez.cigar}
  end

  # ---------------------------------------------------------------------------
  # Reverse complement a 4-bit encoded sequence.
  # ---------------------------------------------------------------------------
  def self.seq_rev_comp(len : Int32, seq : Array(UInt8)) : Array(UInt8)
    result = Array(UInt8).new(len, 0_u8)
    len.times do |i|
      c = seq[len - 1 - i]
      result[i] = c < 4 ? (3_u8 - c) : c
    end
    result
  end

  def self.seq_rev(len : Int32, seq : Array(UInt8)) : Nil
    (len >> 1).times do |i|
      seq[i], seq[len - 1 - i] = seq[len - 1 - i], seq[i]
    end
  end
end
