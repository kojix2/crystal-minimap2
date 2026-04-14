module Minimap2
  # ---------------------------------------------------------------------------
  # Fundamental pair type: (hash<<8|span, rid<<32|pos<<1|strand)
  # Mirrors C's  mm128_t  { uint64_t x, y; }
  # ---------------------------------------------------------------------------
  record Mm128, x : UInt64, y : UInt64 do
    UINT64_MAX = 0xffff_ffff_ffff_ffff_u64

    def self.max
      new(UINT64_MAX, UINT64_MAX)
    end
  end

  # ---------------------------------------------------------------------------
  # One reference sequence inside the index
  # ---------------------------------------------------------------------------
  class IdxSeq
    property name : String
    property offset : UInt64 # offset in the packed-sequence array
    property len : UInt32
    property is_alt : UInt32

    def initialize(@name, @offset, @len, @is_alt = 0_u32)
    end
  end

  # ---------------------------------------------------------------------------
  # Extra alignment information (CIGAR + DP scores) attached to MmReg1
  # ---------------------------------------------------------------------------
  class MmExtra
    property capacity : UInt32
    property dp_score : Int32
    property dp_max : Int32
    property dp_max2 : Int32
    property dp_max0 : Int32
    property n_ambi : UInt32       # number of ambiguous bases
    property trans_strand : UInt32 # transcript strand: 0=unknown, 1=+, 2=-
    property cigar : Array(UInt32) # each element: length<<4 | op

    def initialize
      @capacity = 0_u32
      @dp_score = @dp_max = @dp_max2 = @dp_max0 = 0
      @n_ambi = @trans_strand = 0_u32
      @cigar = [] of UInt32
    end

    def n_cigar : Int32
      @cigar.size
    end
  end

  # ---------------------------------------------------------------------------
  # One alignment hit  (mirrors mm_reg1_t)
  # ---------------------------------------------------------------------------
  class MmReg1
    property id : Int32
    property cnt : Int32    # number of minimizers
    property rid : Int32    # reference index
    property score : Int32  # DP alignment score
    property qs : Int32     # query start
    property qe : Int32     # query end
    property rs : Int32     # reference start
    property re : Int32     # reference end
    property parent : Int32 # parent id (self if primary)
    property subsc : Int32  # best alternate mapping score
    property a_off : Int32  # offset in the a[] array (internal)
    property mlen : Int32   # matching bases
    property blen : Int32   # alignment block length
    property n_sub : Int32  # number of suboptimal mappings
    property score0 : Int32 # initial chaining score

    # Packed bitfields (see accessor methods below)
    property mapq : UInt32
    property split : UInt32
    property rev : Bool
    property inv : Bool
    property sam_pri : Bool
    property proper_frag : Bool
    property pe_thru : Bool
    property seg_split : Bool
    property seg_id : UInt32
    property split_inv : Bool
    property is_alt : Bool
    property strand_retained : Bool
    property is_spliced : Bool

    property hash : UInt32
    property div : Float32
    property p : MmExtra? # nil when CIGAR not requested

    def initialize
      @id = @cnt = @rid = @score = 0
      @qs = @qe = @rs = @re = 0
      @parent = @subsc = @a_off = 0
      @mlen = @blen = @n_sub = @score0 = 0
      @mapq = @split = @seg_id = 0_u32
      @rev = @inv = @sam_pri = @proper_frag = false
      @pe_thru = @seg_split = @split_inv = @is_alt = false
      @strand_retained = @is_spliced = false
      @hash = 0_u32
      @div = 0.0_f32
      @p = nil
    end

    def strand : Char
      @rev ? '-' : '+'
    end
  end

  # ---------------------------------------------------------------------------
  # Index-construction options  (mirrors mm_idxopt_t)
  # ---------------------------------------------------------------------------
  class MmIdxOpt
    property k : Int32    # k-mer size
    property w : Int32    # window size
    property flag : Int32 # MM_I_* flags
    property bucket_bits : Int32
    property mini_batch_size : Int64
    property batch_size : UInt64

    def initialize
      @k = 15
      @w = 10
      @flag = 0
      @bucket_bits = 14
      @mini_batch_size = 50_000_000_i64
      @batch_size = 8_000_000_000_u64
    end
  end

  # ---------------------------------------------------------------------------
  # Mapping / alignment options  (mirrors mm_mapopt_t)
  # ---------------------------------------------------------------------------
  class MmMapOpt
    property flag : Int64
    property seed : Int32
    property sdust_thres : Int32
    property max_qlen : Int32
    property bw : Int32
    property bw_long : Int32
    property max_gap : Int32
    property max_gap_ref : Int32
    property max_frag_len : Int32
    property max_chain_skip : Int32
    property max_chain_iter : Int32
    property min_cnt : Int32
    property min_chain_score : Int32
    property chain_gap_scale : Float32
    property chain_skip_scale : Float32
    property rmq_size_cap : Int32
    property rmq_inner_dist : Int32
    property rmq_rescue_size : Int32
    property rmq_rescue_ratio : Float32
    property mask_level : Float32
    property mask_len : Int32
    property pri_ratio : Float32
    property best_n : Int32
    property alt_drop : Float32
    property a : Int32 # match score
    property b : Int32 # mismatch penalty
    property q : Int32 # gap open
    property e : Int32 # gap extend
    property q2 : Int32
    property e2 : Int32
    property transition : Int32
    property sc_ambi : Int32
    property noncan : Int32
    property junc_bonus : Int32
    property junc_pen : Int32
    property zdrop : Int32
    property zdrop_inv : Int32
    property end_bonus : Int32
    property min_dp_max : Int32
    property min_ksw_len : Int32
    property anchor_ext_len : Int32
    property anchor_ext_shift : Int32
    property max_clip_ratio : Float32
    property rank_min_len : Int32
    property rank_frac : Float32
    property pe_ori : Int32
    property pe_bonus : Int32
    property jump_min_match : Int32
    property mid_occ_frac : Float32
    property q_occ_frac : Float32
    property min_mid_occ : Int32
    property max_mid_occ : Int32
    property mid_occ : Int32
    property max_occ : Int32
    property max_max_occ : Int32
    property occ_dist : Int32
    property mini_batch_size : Int64
    property max_sw_mat : Int64
    property cap_kalloc : Int64
    property split_prefix : String?

    def initialize
      @flag = 0_i64
      @seed = 0
      @sdust_thres = 0
      @max_qlen = 0
      @bw = 0; @bw_long = 0
      @max_gap = 0; @max_gap_ref = 0
      @max_frag_len = 0
      @max_chain_skip = 0; @max_chain_iter = 0
      @min_cnt = 0; @min_chain_score = 0
      @chain_gap_scale = 0.0_f32; @chain_skip_scale = 0.0_f32
      @rmq_size_cap = 0; @rmq_inner_dist = 0; @rmq_rescue_size = 0
      @rmq_rescue_ratio = 0.0_f32
      @mask_level = 0.0_f32; @mask_len = 0
      @pri_ratio = 0.0_f32; @best_n = 0
      @alt_drop = 0.0_f32
      @a = 0; @b = 0; @q = 0; @e = 0; @q2 = 0; @e2 = 0
      @transition = 0; @sc_ambi = 0
      @noncan = 0; @junc_bonus = 0; @junc_pen = 0
      @zdrop = 0; @zdrop_inv = 0; @end_bonus = 0
      @min_dp_max = 0; @min_ksw_len = 0
      @anchor_ext_len = 0; @anchor_ext_shift = 0
      @max_clip_ratio = 0.0_f32
      @rank_min_len = 0; @rank_frac = 0.0_f32
      @pe_ori = 0; @pe_bonus = 0
      @jump_min_match = 0
      @mid_occ_frac = 0.0_f32; @q_occ_frac = 0.0_f32
      @min_mid_occ = 0; @max_mid_occ = 0; @mid_occ = 0
      @max_occ = 0; @max_max_occ = 0; @occ_dist = 0
      @mini_batch_size = 0_i64; @max_sw_mat = 0_i64; @cap_kalloc = 0_i64
      @split_prefix = nil
    end
  end

  # ---------------------------------------------------------------------------
  # Per-thread scratch buffer  (mirrors mm_tbuf_t)
  # rep_len: repetitive region length in query
  # frag_gap: for fragment mode
  # ---------------------------------------------------------------------------
  class MmTbuf
    property rep_len : Int32
    property frag_gap : Int32

    def initialize
      @rep_len = 0
      @frag_gap = 0
    end
  end

  # ---------------------------------------------------------------------------
  # Seed match  (mirrors mm_seed_t, internal to mapping)
  # ---------------------------------------------------------------------------
  struct MmSeed
    property n : Int32       # number of reference hits
    property q_pos : UInt32  # query position (pos<<1|strand)
    property q_span : UInt32 # k-mer span on query
    property flt : Bool      # filtered (too frequent)
    property seg_id : UInt32 # segment id (for paired-end)
    property is_tandem : Bool

    # reference hits: array of (rid<<32|pos) packed into UInt64
    property cr : Array(UInt64)

    def initialize(@n, @q_pos, @q_span, @cr, @seg_id = 0_u32,
                   @flt = false, @is_tandem = false)
    end
  end

  # ---------------------------------------------------------------------------
  # Segment  (mirrors mm_seg_t, internal)
  # ---------------------------------------------------------------------------
  class MmSeg
    property n_u : Int32
    property n_a : Int32
    property u : Array(UInt64)
    property a : Array(Mm128)

    def initialize
      @n_u = @n_a = 0
      @u = [] of UInt64
      @a = [] of Mm128
    end
  end
end
