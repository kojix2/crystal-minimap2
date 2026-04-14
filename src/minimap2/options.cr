module Minimap2
  # ---------------------------------------------------------------------------
  # Option initialisation and preset handling
  # Port of options.c
  # ---------------------------------------------------------------------------

  # Initialise MmIdxOpt with default values (mirrors mm_idxopt_init).
  def self.idxopt_init(opt : MmIdxOpt) : Nil
    opt.k               = 15
    opt.w               = 10
    opt.flag            = 0
    opt.bucket_bits     = 14
    opt.mini_batch_size = 50_000_000_i64
    opt.batch_size      = 8_000_000_000_u64
  end

  # Initialise MmMapOpt with default values (mirrors mm_mapopt_init).
  def self.mapopt_init(opt : MmMapOpt) : Nil
    opt.seed            = 11
    opt.mid_occ_frac    = 2e-4_f32
    opt.min_mid_occ     = 10
    opt.max_mid_occ     = 1_000_000
    opt.sdust_thres     = 0
    opt.q_occ_frac      = 0.01_f32

    opt.min_cnt         = 3
    opt.min_chain_score = 40
    opt.bw              = 500
    opt.bw_long         = 20_000
    opt.max_gap         = 5_000
    opt.max_gap_ref     = -1
    opt.max_chain_skip  = 25
    opt.max_chain_iter  = 5_000
    opt.rmq_inner_dist  = 1_000
    opt.rmq_size_cap    = 100_000
    opt.rmq_rescue_size = 1_000
    opt.rmq_rescue_ratio = 0.1_f32
    opt.chain_gap_scale  = 0.8_f32
    opt.chain_skip_scale = 0.0_f32
    opt.max_max_occ      = 4_095
    opt.occ_dist         = 500

    opt.mask_level  = 0.5_f32
    opt.mask_len    = Int32::MAX
    opt.pri_ratio   = 0.8_f32
    opt.best_n      = 5

    opt.alt_drop    = 0.15_f32

    opt.a  = 2; opt.b = 4; opt.q = 4; opt.e = 2; opt.q2 = 24; opt.e2 = 1
    opt.transition  = 0
    opt.sc_ambi     = 1
    opt.zdrop       = 400; opt.zdrop_inv = 200
    opt.end_bonus   = -1
    opt.min_dp_max  = opt.min_chain_score * opt.a
    opt.min_ksw_len = 200
    opt.anchor_ext_len   = 20; opt.anchor_ext_shift = 6
    opt.max_clip_ratio   = 1.0_f32
    opt.mini_batch_size  = 500_000_000_i64
    opt.max_sw_mat       = 100_000_000_i64
    opt.cap_kalloc       = 500_000_000_i64

    opt.rank_min_len = 500
    opt.rank_frac    = 0.9_f32

    opt.pe_ori   = 0   # FF
    opt.pe_bonus = 33

    opt.jump_min_match = 3
  end

  # Apply *preset* on top of already-initialised *io* and *mo*.
  # Mirrors mm_set_opt().  When *preset* is nil the options are first
  # reset to defaults (equivalent to calling mm_set_opt(NULL,...) in C).
  # Returns 0 on success, -1 if the preset is unknown.
  def self.set_opt(preset : String?, io : MmIdxOpt, mo : MmMapOpt) : Int32
    if preset.nil?
      idxopt_init(io)
      mapopt_init(mo)
      return 0
    end

    case preset
    when "lr", "map-ont"
      # same as default — nothing to change
    when "ava-ont"
      io.flag = 0; io.k = 15; io.w = 5
      mo.flag |= F_ALL_CHAINS | F_NO_DIAG | F_NO_DUAL | F_NO_LJOIN
      mo.min_chain_score = 100; mo.pri_ratio = 0.0_f32; mo.max_chain_skip = 25
      mo.bw = mo.bw_long = 2_000
      mo.occ_dist = 0
    when "map10k", "map-pb"
      io.flag |= I_HPC; io.k = 19
    when "ava-pb"
      io.flag |= I_HPC; io.k = 19; io.w = 5
      mo.flag |= F_ALL_CHAINS | F_NO_DIAG | F_NO_DUAL | F_NO_LJOIN
      mo.min_chain_score = 100; mo.pri_ratio = 0.0_f32; mo.max_chain_skip = 25
      mo.bw_long = mo.bw
      mo.occ_dist = 0
    when "lr:hq", "map-hifi", "map-ccs"
      io.flag = 0; io.k = 19; io.w = 19
      mo.max_gap = 10_000
      mo.min_mid_occ = 50; mo.max_mid_occ = 500
      if preset == "map-hifi" || preset == "map-ccs"
        mo.a = 1; mo.b = 4; mo.q = 6; mo.q2 = 26; mo.e = 2; mo.e2 = 1
        mo.min_dp_max = 200
      end
    when "lr:hqae"
      io.flag = 0; io.k = 25; io.w = 51
      mo.flag |= F_RMQ
      mo.min_mid_occ = 50; mo.max_mid_occ = 500
      mo.rmq_inner_dist = 5_000
      mo.occ_dist = 200
      mo.best_n = 100
      mo.chain_gap_scale = 5.0_f32
    when "map-iclr-prerender"
      io.flag = 0; io.k = 15
      mo.b = 6; mo.transition = 1
      mo.q = 10; mo.q2 = 50
    when "map-iclr"
      io.flag = 0; io.k = 19
      mo.b = 6; mo.transition = 4
      mo.q = 10; mo.q2 = 50
    when .starts_with?("asm")
      io.flag = 0; io.k = 19; io.w = 19
      mo.bw = 1_000; mo.bw_long = 100_000
      mo.max_gap = 10_000
      mo.flag |= F_RMQ
      mo.min_mid_occ = 50; mo.max_mid_occ = 500
      mo.min_dp_max = 200
      mo.best_n = 50
      case preset
      when "asm5"
        mo.a = 1; mo.b = 19; mo.q = 39; mo.q2 = 81; mo.e = 3; mo.e2 = 1
        mo.zdrop = mo.zdrop_inv = 200
      when "asm10"
        mo.a = 1; mo.b = 9; mo.q = 16; mo.q2 = 41; mo.e = 2; mo.e2 = 1
        mo.zdrop = mo.zdrop_inv = 200
      when "asm20"
        mo.a = 1; mo.b = 4; mo.q = 6; mo.q2 = 26; mo.e = 2; mo.e2 = 1
        mo.zdrop = mo.zdrop_inv = 200
        io.w = 10
      else
        return -1
      end
    when "short", "sr"
      io.flag = 0; io.k = 21; io.w = 11
      mo.flag |= F_SR | F_FRAG_MODE | F_NO_PRINT_2ND | F_2_IO_THREADS | F_HEAP_SORT
      mo.pe_ori = 0 << 1 | 1  # FR
      mo.a = 2; mo.b = 8; mo.q = 12; mo.e = 2; mo.q2 = 24; mo.e2 = 1
      mo.zdrop = mo.zdrop_inv = 100
      mo.end_bonus = 10
      mo.max_frag_len = 800
      mo.max_gap = 100
      mo.bw = mo.bw_long = 100
      mo.pri_ratio = 0.5_f32
      mo.min_cnt = 2
      mo.min_chain_score = 25
      mo.min_dp_max = 40
      mo.best_n = 20
      mo.mid_occ = 1_000
      mo.max_occ = 5_000
      mo.mini_batch_size = 50_000_000_i64
    when "splice", "splice:hq", "splice:sr", "cdna"
      io.flag = 0; io.k = 15; io.w = 5
      mo.flag |= F_SPLICE | F_SPLICE_FOR | F_SPLICE_REV | F_SPLICE_FLANK
      mo.max_sw_mat = 0_i64
      mo.max_gap = 2_000; mo.max_gap_ref = mo.bw = mo.bw_long = 200_000
      mo.a = 1; mo.b = 2; mo.q = 2; mo.e = 1; mo.q2 = 32; mo.e2 = 0
      mo.noncan = 9
      mo.junc_bonus = 9
      mo.junc_pen = 5
      mo.zdrop = 200; mo.zdrop_inv = 100
      case preset
      when "splice:hq"
        mo.noncan = 5; mo.b = 4; mo.q = 6; mo.q2 = 24
      when "splice:sr"
        mo.flag |= F_NO_PRINT_2ND | F_2_IO_THREADS | F_HEAP_SORT | F_FRAG_MODE | F_WEAK_PAIRING | F_SR_RNA
        mo.noncan = 5; mo.b = 4; mo.q = 6; mo.q2 = 24
        mo.min_chain_score = 25
        mo.min_dp_max = 40
        mo.min_ksw_len = 20
        mo.pe_ori = 0 << 1 | 1  # FR
        mo.best_n = 10
        mo.mini_batch_size = 100_000_000_i64
      end
    else
      return -1
    end
    0
  end

  # Update mid_occ in *opt* based on the index *mi* (mirrors mm_mapopt_update).
  def self.mapopt_update(opt : MmMapOpt, mi : MmIdx) : Nil
    if (opt.flag & F_SPLICE_FOR) != 0 || (opt.flag & F_SPLICE_REV) != 0
      opt.flag |= F_SPLICE
    end
    if opt.mid_occ <= 0
      opt.mid_occ = mi.cal_max_occ(opt.mid_occ_frac)
      opt.mid_occ = opt.min_mid_occ if opt.mid_occ < opt.min_mid_occ
      if opt.max_mid_occ > opt.min_mid_occ && opt.mid_occ > opt.max_mid_occ
        opt.mid_occ = opt.max_mid_occ
      end
    end
    opt.bw_long = opt.bw if opt.bw_long < opt.bw
    if Minimap2.verbose >= 3
      t = Minimap2.realtime - Minimap2.realtime0
      STDERR.printf("[M::mapopt_update::%.3f*%.2f] mid_occ = %d\n", t, Minimap2.cputime / t, opt.mid_occ)
    end
  end

  # Set bw/bw_long/max_gap_ref based on intron length (mirrors mm_mapopt_max_intron_len).
  def self.mapopt_max_intron_len(opt : MmMapOpt, max_intron_len : Int32) : Nil
    if (opt.flag & F_SPLICE) != 0 && max_intron_len > 0
      opt.max_gap_ref = opt.bw = opt.bw_long = max_intron_len
    end
  end

  # Maximum supplementary-chain bonus (mirrors mm_max_spsc_bonus).
  def self.max_spsc_bonus(mo : MmMapOpt) : Int32
    max_sc = (mo.q2 + 1) / 2 - 1
    max_sc = mo.q2 - mo.q if mo.q2 - mo.q > max_sc
    max_sc
  end

  # Validate option combination; returns 0 on success (mirrors mm_check_opt).
  def self.check_opt(io : MmIdxOpt, mo : MmMapOpt) : Int32
    if mo.bw > mo.bw_long
      STDERR.puts "[ERROR] with '-rNUM1,NUM2', NUM1 (#{mo.bw}) can't be larger than NUM2 (#{mo.bw_long})" if Minimap2.verbose >= 1
      return -8
    end
    if (mo.flag & F_RMQ) != 0 && (mo.flag & (F_SR | F_SPLICE)) != 0
      STDERR.puts "[ERROR] --rmq doesn't work with --sr or --splice" if Minimap2.verbose >= 1
      return -7
    end
    if mo.split_prefix && (mo.flag & (F_OUT_CS | F_OUT_MD)) != 0
      STDERR.puts "[ERROR] --cs or --MD doesn't work with --split-prefix" if Minimap2.verbose >= 1
      return -6
    end
    if io.k <= 0 || io.w <= 0
      STDERR.puts "[ERROR] -k and -w must be positive" if Minimap2.verbose >= 1
      return -5
    end
    if mo.best_n < 0
      STDERR.puts "[ERROR] -N must be no less than 0" if Minimap2.verbose >= 1
      return -4
    end
    if mo.best_n == 0 && Minimap2.verbose >= 2
      STDERR.puts "[WARNING] '-N 0' reduces mapping accuracy. Please use '--secondary=no' instead."
    end
    if mo.pri_ratio < 0.0_f32 || mo.pri_ratio > 1.0_f32
      STDERR.puts "[ERROR] -p must be within 0 and 1" if Minimap2.verbose >= 1
      return -4
    end
    if (mo.flag & F_FOR_ONLY) != 0 && (mo.flag & F_REV_ONLY) != 0
      STDERR.puts "[ERROR] --for-only and --rev-only can't be applied at the same time" if Minimap2.verbose >= 1
      return -3
    end
    if mo.e <= 0 || mo.q <= 0
      STDERR.puts "[ERROR] -O and -E must be positive" if Minimap2.verbose >= 1
      return -1
    end
    if (mo.q != mo.q2 || mo.e != mo.e2) && !(mo.e > mo.e2 && mo.q + mo.e < mo.q2 + mo.e2)
      STDERR.puts "[ERROR] dual gap penalties violating E1>E2 and O1+E1<O2+E2" if Minimap2.verbose >= 1
      return -2
    end
    if (mo.q + mo.e) + (mo.q2 + mo.e2) > 127
      STDERR.puts "[ERROR] scoring system violating ({-O}+{-E})+({-O2}+{-E2}) <= 127" if Minimap2.verbose >= 1
      return -1
    end
    if mo.sc_ambi < 0 || mo.sc_ambi >= mo.b
      STDERR.puts "[ERROR] --score-N should be within [0,{-B})" if Minimap2.verbose >= 1
      return -1
    end
    if mo.zdrop < mo.zdrop_inv
      STDERR.puts "[ERROR] Z-drop should not be less than inversion-Z-drop" if Minimap2.verbose >= 1
      return -5
    end
    if (mo.flag & F_NO_PRINT_2ND) != 0 && (mo.flag & F_ALL_CHAINS) != 0
      STDERR.puts "[ERROR] -X/-P and --secondary=no can't be applied at the same time" if Minimap2.verbose >= 1
      return -5
    end
    0
  end
end
