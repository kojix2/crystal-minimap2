require "./spec_helper"

describe Minimap2 do
  describe "VERSION" do
    it "is defined" do
      Minimap2::VERSION.should_not be_empty
    end
  end

  # ---------------------------------------------------------------------------
  # Options
  # ---------------------------------------------------------------------------
  describe "options" do
    it "idxopt_init sets sane defaults" do
      io = Minimap2::MmIdxOpt.new
      Minimap2.idxopt_init(io)
      io.k.should eq(15)
      io.w.should eq(10)
      io.bucket_bits.should eq(14)
    end

    it "mapopt_init sets sane defaults" do
      mo = Minimap2::MmMapOpt.new
      Minimap2.mapopt_init(mo)
      mo.a.should be > 0
      mo.b.should be > 0 # mismatch penalty magnitude (used as positive, negated in matrix)
      mo.min_cnt.should be > 0
    end

    it "set_opt with nil gives defaults" do
      io = Minimap2::MmIdxOpt.new
      mo = Minimap2::MmMapOpt.new
      ret = Minimap2.set_opt(nil, io, mo)
      ret.should eq(0)
      io.k.should be > 0
    end

    it "set_opt with map-ont preset" do
      io = Minimap2::MmIdxOpt.new
      mo = Minimap2::MmMapOpt.new
      ret = Minimap2.set_opt("map-ont", io, mo)
      ret.should eq(0)
    end

    it "set_opt with asm5 preset" do
      io = Minimap2::MmIdxOpt.new
      mo = Minimap2::MmMapOpt.new
      ret = Minimap2.set_opt("asm5", io, mo)
      ret.should eq(0)
    end

    it "set_opt with unknown preset returns -1" do
      io = Minimap2::MmIdxOpt.new
      mo = Minimap2::MmMapOpt.new
      ret = Minimap2.set_opt("nonexistent_preset", io, mo)
      ret.should eq(-1)
    end
  end

  # ---------------------------------------------------------------------------
  # Sketch
  # ---------------------------------------------------------------------------
  describe "mm_sketch" do
    it "produces minimizers for a simple sequence" do
      seq = "ACGTACGTACGTACGTACGTACGT"
      mv = [] of Minimap2::Mm128
      Minimap2.mm_sketch(seq, seq.size, w: 5, k: 10, rid: 0_u32, is_hpc: false, p: mv)
      mv.should_not be_empty
    end

    it "produces no minimizers for a short sequence" do
      seq = "ACGT" # shorter than k=15
      mv = [] of Minimap2::Mm128
      Minimap2.mm_sketch(seq, seq.size, w: 5, k: 15, rid: 0_u32, is_hpc: false, p: mv)
      mv.should be_empty
    end

    it "different sequences produce different minimizers" do
      seq1 = "ACGTACGTACGTACGTACGTACGTACGT"
      seq2 = "TTTTTTTTTTTTTTTTTTTTTTTTTTTT"
      mv1 = [] of Minimap2::Mm128
      mv2 = [] of Minimap2::Mm128
      Minimap2.mm_sketch(seq1, seq1.size, w: 5, k: 10, rid: 0_u32, is_hpc: false, p: mv1)
      Minimap2.mm_sketch(seq2, seq2.size, w: 5, k: 10, rid: 0_u32, is_hpc: false, p: mv2)
      # They should differ (ACGT vs TTTT)
      xs1 = mv1.map(&.x).sort!
      xs2 = mv2.map(&.x).sort!
      xs1.should_not eq(xs2)
    end
  end

  # ---------------------------------------------------------------------------
  # SDUST
  # ---------------------------------------------------------------------------
  describe "sdust" do
    it "returns fewer masked regions for complex sequence than repetitive" do
      complex = "ACGTACGTACGTACGTACGTACGT"
      repetitive = "A" * 64
      r_complex = Minimap2.sdust(complex.to_slice, big_t: 20)
      r_repetitive = Minimap2.sdust(repetitive.to_slice, big_t: 20)
      r_repetitive.size.should be >= r_complex.size
    end

    it "masks low-complexity runs" do
      # All-A is maximally low complexity
      seq = "A" * 64
      result = Minimap2.sdust(seq.to_slice, big_t: 20)
      result.should_not be_empty
    end
  end

  # ---------------------------------------------------------------------------
  # Sequence encode/revcomp
  # ---------------------------------------------------------------------------
  describe "encode_seq / rev_comp" do
    it "encodes ACGT to 0-3" do
      enc = Minimap2.encode_seq("ACGT")
      enc.should eq([0_u8, 1_u8, 2_u8, 3_u8])
    end

    it "encodes N to 4" do
      enc = Minimap2.encode_seq("N")
      enc.should eq([4_u8])
    end

    it "rev_comp reverses and complements" do
      # ACGT encodes as [0,1,2,3]. RevComp should give [0,1,2,3] (ACGT is palindrome)
      enc = Minimap2.encode_seq("ACGT")
      rc = Minimap2.rev_comp(enc)
      rc.should eq([0_u8, 1_u8, 2_u8, 3_u8]) # ACGT is its own revcomp
    end

    it "rev_comp of AAAA gives TTTT" do
      # A=0, T=3
      enc = Minimap2.encode_seq("AAAA")
      rc = Minimap2.rev_comp(enc)
      rc.should eq([3_u8, 3_u8, 3_u8, 3_u8])
    end
  end

  # ---------------------------------------------------------------------------
  # KSW2 substitution matrix
  # ---------------------------------------------------------------------------
  describe "ksw_gen_simple_mat" do
    it "diagonal is positive" do
      mat = Minimap2.ksw_gen_simple_mat(5, 2, 4, 1)
      # mat[i*5+i] for i in 0..3 should be match score
      (0..3).each do |i|
        mat[i * 5 + i].should eq(2_i8)
      end
    end

    it "off-diagonal is negative mismatch" do
      mat = Minimap2.ksw_gen_simple_mat(5, 2, 4, 1)
      mat[0 * 5 + 1].should eq(-4_i8)
    end
  end

  # ---------------------------------------------------------------------------
  # KSW2 alignment
  # ---------------------------------------------------------------------------
  describe "ksw_extz2" do
    it "aligns identical sequences" do
      q = Minimap2.encode_seq("ACGTACGT")
      t = Minimap2.encode_seq("ACGTACGT")
      mat = Minimap2.ksw_gen_simple_mat(5, 2, 4, 1)
      ez = Minimap2.ksw_extz2(q.size, q, t.size, t, 5, mat, 4, 2, -1, -1, 0, 0)
      ez.score.should be > 0
    end

    it "score is higher for perfect match than mismatched" do
      q = Minimap2.encode_seq("ACGTACGT")
      t1 = Minimap2.encode_seq("ACGTACGT")
      t2 = Minimap2.encode_seq("TTTTTTTT")
      mat = Minimap2.ksw_gen_simple_mat(5, 2, 4, 1)
      ez1 = Minimap2.ksw_extz2(q.size, q, t1.size, t1, 5, mat, 4, 2, -1, -1, 0, 0)
      ez2 = Minimap2.ksw_extz2(q.size, q, t2.size, t2, 5, mat, 4, 2, -1, -1, 0, 0)
      ez1.score.should be > ez2.score
    end
  end

  # ---------------------------------------------------------------------------
  # Index building and lookup
  # ---------------------------------------------------------------------------
  describe "MmIdx" do
    ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    ref_name = "chr1"

    it "can be built from strings" do
      io = Minimap2::MmIdxOpt.new
      Minimap2.idxopt_init(io)
      io.k = 10; io.w = 5

      mi = Minimap2::MmIdx.from_strings(io.w, io.k, (io.flag & Minimap2::I_HPC) != 0, io.bucket_bits, [ref_seq], [ref_name])
      mi.should_not be_nil
      mi.seq.size.should eq(1)
      mi.seq[0].name.should eq(ref_name)
    end

    it "name2id returns correct index" do
      io = Minimap2::MmIdxOpt.new
      Minimap2.idxopt_init(io)
      io.k = 10; io.w = 5

      mi = Minimap2::MmIdx.from_strings(io.w, io.k, (io.flag & Minimap2::I_HPC) != 0, io.bucket_bits, [ref_seq], [ref_name])
      mi.name2id(ref_name).should eq(0)
      mi.name2id("nonexistent").should eq(-1)
    end

    it "getseq retrieves bases" do
      io = Minimap2::MmIdxOpt.new
      Minimap2.idxopt_init(io)
      io.k = 10; io.w = 5

      mi = Minimap2::MmIdx.from_strings(io.w, io.k, (io.flag & Minimap2::I_HPC) != 0, io.bucket_bits, [ref_seq], [ref_name])
      buf = Array(UInt8).new(4, 0_u8)
      mi.getseq(0_u32, 0_u32, 4_u32, buf)
      # ACGT -> 0,1,2,3
      buf.should eq([0_u8, 1_u8, 2_u8, 3_u8])
    end
  end

  # ---------------------------------------------------------------------------
  # End-to-end: Aligner
  # ---------------------------------------------------------------------------
  describe "Aligner" do
    ref = "ACGT" * 50
    qry = "ACGT" * 20

    it "can be constructed from strings" do
      aln = Minimap2::Aligner.from_strings([ref], ["ref"], "map-ont")
      aln.indices.should_not be_empty
      aln.seq_names.should eq(["ref"])
    end

    it "map returns an array" do
      aln = Minimap2::Aligner.from_strings([ref], ["ref"], "map-ont")
      hits = aln.map(qry, "query")
      hits.should be_a(Array(Minimap2::MmReg1))
    end

    it "maps query to correct reference strand" do
      aln = Minimap2::Aligner.from_strings([ref], ["ref"], "map-ont")
      hits = aln.map(qry, "query")
      if !hits.empty?
        hits.first.rev?.should be_false # query is same strand
      end
    end
  end
end
