module Minimap2
  # ---------------------------------------------------------------------------
  # Library version
  # ---------------------------------------------------------------------------
  LIB_VERSION = "2.30"

  # ---------------------------------------------------------------------------
  # Mapping / output flags  (MM_F_*)
  # ---------------------------------------------------------------------------
  F_NO_DIAG       =         0x001_i64
  F_NO_DUAL       =         0x002_i64
  F_CIGAR         =         0x004_i64
  F_OUT_SAM       =         0x008_i64
  F_NO_QUAL       =         0x010_i64
  F_OUT_CG        =         0x020_i64
  F_OUT_CS        =         0x040_i64
  F_SPLICE        =         0x080_i64
  F_SPLICE_FOR    =         0x100_i64
  F_SPLICE_REV    =         0x200_i64
  F_NO_LJOIN      =         0x400_i64
  F_OUT_CS_LONG   =         0x800_i64
  F_SR            =        0x1000_i64
  F_FRAG_MODE     =        0x2000_i64
  F_NO_PRINT_2ND  =        0x4000_i64
  F_2_IO_THREADS  =        0x8000_i64
  F_LONG_CIGAR    =       0x10000_i64
  F_INDEPEND_SEG  =       0x20000_i64
  F_SPLICE_FLANK  =       0x40000_i64
  F_SOFTCLIP      =       0x80000_i64
  F_FOR_ONLY      =      0x100000_i64
  F_REV_ONLY      =      0x200000_i64
  F_HEAP_SORT     =      0x400000_i64
  F_ALL_CHAINS    =      0x800000_i64
  F_OUT_MD        =     0x1000000_i64
  F_COPY_COMMENT  =     0x2000000_i64
  F_EQX           =     0x4000000_i64
  F_PAF_NO_HIT    =     0x8000000_i64
  F_NO_END_FLT    =    0x10000000_i64
  F_HARD_MLEVEL   =    0x20000000_i64
  F_SAM_HIT_ONLY  =    0x40000000_i64
  F_RMQ           =    0x80000000_i64
  F_QSTRAND       =   0x100000000_i64
  F_NO_INV        =   0x200000000_i64
  F_NO_HASH_NAME  =   0x400000000_i64
  F_SPLICE_OLD    =   0x800000000_i64
  F_SECONDARY_SEQ =  0x1000000000_i64
  F_OUT_DS        =  0x2000000000_i64
  F_WEAK_PAIRING  =  0x4000000000_i64
  F_SR_RNA        =  0x8000000000_i64
  F_OUT_JUNC      = 0x10000000000_i64

  # ---------------------------------------------------------------------------
  # Index flags  (MM_I_*)
  # ---------------------------------------------------------------------------
  I_HPC     = 0x1
  I_NO_SEQ  = 0x2
  I_NO_NAME = 0x4

  # ---------------------------------------------------------------------------
  # CIGAR operation codes
  # ---------------------------------------------------------------------------
  CIGAR_MATCH      = 0
  CIGAR_INS        = 1
  CIGAR_DEL        = 2
  CIGAR_N_SKIP     = 3
  CIGAR_SOFTCLIP   = 4
  CIGAR_HARDCLIP   = 5
  CIGAR_PADDING    = 6
  CIGAR_EQ_MATCH   = 7
  CIGAR_X_MISMATCH = 8

  CIGAR_STR = "MIDNSHP=XB"

  # ---------------------------------------------------------------------------
  # Internal seed flags
  # ---------------------------------------------------------------------------
  SEED_LONG_JOIN = 1_u64 << 40
  SEED_IGNORE    = 1_u64 << 41
  SEED_TANDEM    = 1_u64 << 42
  SEED_SELF      = 1_u64 << 43

  SEED_SEG_SHIFT = 48
  SEED_SEG_MASK  = 0xff_u64 << 48

  # ---------------------------------------------------------------------------
  # Debug flags
  # ---------------------------------------------------------------------------
  DBG_NO_KALLOC     =  0x1
  DBG_PRINT_QNAME   =  0x2
  DBG_PRINT_SEED    =  0x4
  DBG_PRINT_ALN_SEQ =  0x8
  DBG_PRINT_CHAIN   = 0x10
  DBG_SEED_FREQ     = 0x20

  # ---------------------------------------------------------------------------
  # Index magic
  # ---------------------------------------------------------------------------
  IDX_MAGIC = "MMI\x02"
  MAX_SEG   = 255

  # ---------------------------------------------------------------------------
  # Internal parent sentinels
  # ---------------------------------------------------------------------------
  PARENT_UNSET   = -1
  PARENT_TMP_PRI = -2

  # ---------------------------------------------------------------------------
  # Junction flags
  # ---------------------------------------------------------------------------
  JUNC_ANNO = 0x1
  JUNC_MISC = 0x2
end
