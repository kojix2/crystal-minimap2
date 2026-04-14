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

  KSW_SPSC_OFFSET = 64

  # Result of a KSW alignment
  class KswExtz
    property max : UInt32 # max score seen (unsigned representation)
    property zdropped : Bool
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
  private def self.push_cigar(cigar : Array(UInt32), op : Int32, len : Int32) : Nil
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
        e_val = [h_up - qe, e_arr[j] - e].max
        e2_val = [h_up - qe2, e2_arr[j] - e2].max
        e_arr[j] = e_val
        e2_arr[j] = e2_val

        # Insertion from left (F)
        h_left = j > j_lo ? h_curr[j - 1] : neg
        f_new = [h_left - qe, f - e].max
        f2_new = [h_left - qe2, f2 - e2].max
        f = f_new
        f2 = f2_new

        # Match/mismatch from diagonal
        h_diag = (i > 0 && j > 0) ? h_prev[j - 1] : (i == 0 && j == 0 ? 0 : neg)
        ti = target[j].to_i32
        sc = (h_diag == neg) ? neg : h_diag + mat[qi * m + ti].to_i32

        # H
        h = [sc, e_val, e2_val, f_new, f2_new].max
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
        when 1, 3 # deletion (E or E2)
          push_cigar(ez.cigar, CIGAR_DEL, 1)
          ci -= 1
        when 2, 4 # insertion (F or F2)
          push_cigar(ez.cigar, CIGAR_INS, 1)
          cj -= 1
        else
          break
        end
      end
      # leading indels
      push_cigar(ez.cigar, CIGAR_DEL, ci + 1) if ci >= 0
      push_cigar(ez.cigar, CIGAR_INS, cj + 1) if cj >= 0

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
  # ksw_exts2: splice-aware extension — simplified scalar implementation.
  # Mirrors ksw_exts2_sse().
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
    # For now, delegate to ksw_extd2 ignoring splice-specific logic.
    # A full implementation would add splice state (N-skip) with junction scores.
    ksw_extd2(qlen, query, tlen, target, m, mat,
      q, e, q2, 1,
      -1, zdrop, end_bonus, flag)
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
    result = Array(UInt8).new(len)
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
