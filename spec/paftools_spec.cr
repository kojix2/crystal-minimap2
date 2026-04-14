require "spec"

# Run the pre-built paftools binary and capture stdout + stderr.
PAFTOOLS_BIN = File.expand_path("../bin/paftools", __DIR__)
FIXTURES     = File.expand_path("fixtures", __DIR__)

private def run_paftools(*args) : {Int32, String, String}
  out_io = IO::Memory.new
  err_io = IO::Memory.new
  status = Process.run(PAFTOOLS_BIN, args: args.to_a,
    output: out_io, error: err_io)
  {status.exit_code, out_io.to_s, err_io.to_s}
end

describe "paftools binary" do
  it "exists" do
    File.exists?(PAFTOOLS_BIN).should be_true
  end

  # ── version ──────────────────────────────────────────────────────────────
  describe "version" do
    it "prints version string and exits 0" do
      code, out, _ = run_paftools("version")
      code.should eq(0)
      out.strip.should match(/2\.30/)
    end
  end

  # ── stat ─────────────────────────────────────────────────────────────────
  describe "stat" do
    # Fixture: 2 primary reads + 1 secondary (s2:i: tag)
    # read1: 1000 bp, CIGAR 500M10D490M, NM=10  → 0 substitutions, 1 deletion in [0,50)
    # read2:  800 bp, CIGAR 800M,         NM=5  → 5 substitutions
    # read1 secondary: has s2:i: → counted in n_2nd, skipped from stats

    it "counts sequences, primary and secondary alignments" do
      code, out, _ = run_paftools("stat", "#{FIXTURES}/stat.paf")
      code.should eq(0)
      out.should match(/Number of mapped sequences: 2/)
      out.should match(/Number of primary alignments: 2/)
      out.should match(/Number of secondary alignments: 1/)
    end

    it "counts substitutions correctly" do
      _, out, _ = run_paftools("stat", "#{FIXTURES}/stat.paf")
      out.should match(/Number of substitutions: 5/)
    end

    it "reports zero 65k-CIGAR alignments" do
      _, out, _ = run_paftools("stat", "#{FIXTURES}/stat.paf")
      out.should match(/Number of primary alignments with >65535 CIGAR operations: 0/)
    end

    it "counts deletions by size bucket" do
      _, out, _ = run_paftools("stat", "#{FIXTURES}/stat.paf")
      # read1 has one 10-bp deletion, which falls in [0,50)
      out.should match(/Number of deletions in \[0,50\): 1/)
      out.should match(/Number of deletions in \[50,100\): 0/)
    end

    it "reports mapped bases" do
      _, out, _ = run_paftools("stat", "#{FIXTURES}/stat.paf")
      # read1 qlen=1000, read2 qlen=800 → total 1800; both fully mapped → cov 1800
      out.should match(/Number of bases in mapped sequences: 1800/)
      out.should match(/Number of mapped bases: 1800/)
    end
  end

  # ── sam2paf ───────────────────────────────────────────────────────────────
  describe "sam2paf" do
    # Fixture: @SQ chr1 LN:5000 + read1 flag=0 chr1:101 CIGAR=10M NM=0
    # Expected PAF: read1 10 0 10 + chr1 5000 100 110 10 10 60 tp:A:P NM:i:0 ...

    it "converts a simple SAM to PAF" do
      code, out, _ = run_paftools("sam2paf", "#{FIXTURES}/sam2paf.sam")
      code.should eq(0)
      fields = out.strip.split('\t')
      fields[0].should eq("read1") # qname
      fields[1].should eq("10")    # qlen
      fields[2].should eq("0")     # qs
      fields[3].should eq("10")    # qe
      fields[4].should eq("+")     # strand
      fields[5].should eq("chr1")  # tname
      fields[6].should eq("5000")  # tlen
      fields[7].should eq("100")   # ts  (1-based 101 → 0-based 100)
      fields[8].should eq("110")   # te
      fields[9].should eq("10")    # matches
      fields[10].should eq("10")   # block_len
      fields[11].should eq("60")   # mapq
    end

    it "includes tp, NM and cg tags" do
      _, out, _ = run_paftools("sam2paf", "#{FIXTURES}/sam2paf.sam")
      out.should match(/tp:A:P/)
      out.should match(/NM:i:0/)
      out.should match(/cg:Z:10M/)
    end
  end

  # ── delta2paf ────────────────────────────────────────────────────────────
  describe "delta2paf" do
    # Fixture: a perfect 80-bp alignment of qry1[0,80) onto ref1[0,80)

    it "converts a MUMmer delta to PAF" do
      code, out, _ = run_paftools("delta2paf", "#{FIXTURES}/test.delta")
      code.should eq(0)
      fields = out.strip.split('\t')
      fields[0].should eq("qry1") # qname
      fields[1].should eq("80")   # qlen
      fields[2].should eq("0")    # qs
      fields[3].should eq("80")   # qe
      fields[4].should eq("+")    # strand
      fields[5].should eq("ref1") # tname
      fields[6].should eq("100")  # tlen
      fields[7].should eq("0")    # ts
      fields[8].should eq("80")   # te
      fields[9].should eq("80")   # matches (blen - nm = 80 - 0)
      fields[10].should eq("80")  # block_len
    end

    it "emits NM:i:0 and cg:Z:80M tags" do
      _, out, _ = run_paftools("delta2paf", "#{FIXTURES}/test.delta")
      out.should match(/NM:i:0/)
      out.should match(/cg:Z:80M/)
    end
  end

  # ── longcs2seq ───────────────────────────────────────────────────────────
  describe "longcs2seq" do
    # Fixture cs: =ACGT*ac=ACGTACGT
    # ref sequence:   ACGT + A (ref base of *ac) + ACGTACGT = ACGTAACGTACGT
    # query sequence: ACGT + C (qry base of *ac) + ACGTACGT = ACGTCACGTACGT

    it "reconstructs the reference sequence" do
      code, out, _ = run_paftools("longcs2seq", "#{FIXTURES}/longcs.paf")
      code.should eq(0)
      lines = out.strip.split('\n')
      lines[0].should eq(">ref1_0_13")
      lines[1].should eq("ACGTAACGTACGT")
    end

    it "reconstructs the query sequence with -q" do
      code, out, _ = run_paftools("longcs2seq", "-q", "#{FIXTURES}/longcs.paf")
      code.should eq(0)
      lines = out.strip.split('\n')
      lines[0].should eq(">read1_5_18")
      lines[1].should eq("ACGTCACGTACGT")
    end
  end

  # ── view ─────────────────────────────────────────────────────────────────
  describe "view" do
    # Fixture: read1 vs ref1, cg:Z:20M, AS:i:100

    it "produces lastz-cigar output" do
      code, out, _ = run_paftools("view", "-f", "lastz-cigar", "#{FIXTURES}/view.paf")
      code.should eq(0)
      out.strip.should eq("cigar: read1 0 20 + ref1 0 20 + 100 M 20")
    end
  end

  # ── gff2bed ───────────────────────────────────────────────────────────────
  describe "gff2bed" do
    # Fixture: 2-exon transcript tx1 on chr1+, with a CDS
    # Exons:  [100,200) and [300,400) (0-based, from GTF 1-based 101-200 and 301-400)
    # CDS:    [100,200)

    it "converts GTF exons to BED12" do
      code, out, _ = run_paftools("gff2bed", "#{FIXTURES}/test.gtf")
      code.should eq(0)
      fields = out.strip.split('\t')
      fields[0].should eq("chr1")
      fields[1].should eq("100")                      # BED start (0-based)
      fields[2].should eq("400")                      # BED end
      fields[3].should eq("tx1|protein_coding|geneA") # name
      fields[4].should eq("1000")                     # score
      fields[5].should eq("+")                        # strand
      fields[6].should eq("100")                      # thickStart
      fields[7].should eq("200")                      # thickEnd
      fields[9].should eq("2")                        # block count
      fields[10].should eq("100,100,")                # block sizes
      fields[11].should eq("0,200,")                  # block starts
    end
  end

  # ── bedcov ───────────────────────────────────────────────────────────────
  describe "bedcov" do
    # regions.bed:  chr1:100-200, chr1:300-500
    # target.bed:   chr1:50-250  (overlaps 100 bp with [100,200))
    #               chr1:400-600 (overlaps 100 bp with [300,500))
    # tot_len = 200+200 = 400, hit_len = 100+100 = 200  (50%)

    it "reports coverage fractions on stderr" do
      code, _, err = run_paftools("bedcov",
        "#{FIXTURES}/regions.bed",
        "#{FIXTURES}/target.bed")
      code.should eq(0)
      err.should match(/# target bases: 400/)
      err.should match(/# target bases overlapping regions: 200 \(50\.00%\)/)
    end
  end

  # ── liftover ─────────────────────────────────────────────────────────────
  describe "liftover" do
    # aln.paf:   query1[100,600) ↔ chr1[1000,1500) via 500M (perfect match)
    # query.bed: query1:200-400
    # Lifted:    chr1:1100-1300
    # Use -l 100 to lower the default min_len=50000 filter

    it "lifts BED coordinates through a PAF alignment" do
      code, out, _ = run_paftools("liftover", "-l", "100",
        "#{FIXTURES}/aln.paf",
        "#{FIXTURES}/query.bed")
      code.should eq(0)
      fields = out.strip.split('\t')
      fields[0].should eq("chr1")
      fields[1].should eq("1100")
      fields[2].should eq("1300")
      fields[3].should eq("query1_200_400")
      fields[5].should eq("+")
    end
  end
end
