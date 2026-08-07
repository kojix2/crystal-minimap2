require "spec"

# Run the pre-built paftools binary and capture stdout + stderr.
PAFTOOLS_BIN = ENV["PAFTOOLS_BIN"]? || begin
  exe = {% if flag?(:win32) %}
          "paftools.exe"
        {% else %}
          "paftools"
        {% end %}
  File.expand_path("../bin/#{exe}", __DIR__)
end
FIXTURES = File.expand_path("fixtures", __DIR__)

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
    # Fixture: 2 primary reads + 1 secondary alignment.
    # minimap2 only tags the PRIMARY record with s2:i: (2nd-best chaining
    # score, format.c: `if (r->parent == r->id) ... s2:i:`), so paftools.js
    # treats the *absence* of s2:i: as the marker of a secondary alignment.
    # read1: 1000 bp, CIGAR 500M10D490M, NM=10  → 0 substitutions, 1 deletion in [0,50)
    # read2:  800 bp, CIGAR 800M,         NM=5  → 5 substitutions
    # read1 secondary (tp:A:S, no s2:i:) → counted in n_2nd, skipped from stats

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

  # ── misjoin ──────────────────────────────────────────────────────────────
  describe "misjoin" do
    # q1: chrA→chrB→chrA (inter-chromosomal 'J' x2); q3: +/-/+ bracketed
    # inversion on chrD ('M' x3). Every row carries trailing tp:A:P/cm:i:100
    # tags so the -p full-row test below actually exercises the truncation
    # bug fix. Verified byte-for-byte against k8 paftools.js.

    it "detects inter-chromosomal misjoins and bracketed inversions" do
      code, out, _ = run_paftools("misjoin", "-e", "#{FIXTURES}/misjoin.paf")
      code.should eq(0)
      out.should match(/^J\tq1/m)
      out.should match(/^M\tq3/m)
      out.should match(/# inter-chromosomal misjoins: 2,0/)
      out.should match(/# candidate inversions in the middle: 1,0/)
    end

    it "truncates -e error-report rows to the first 12 PAF columns (tags dropped)" do
      code, out, _ = run_paftools("misjoin", "-e", "#{FIXTURES}/misjoin.paf")
      code.should eq(0)
      j_line = out.lines.find { |l| l.starts_with?("J\tq1") }.not_nil!
      j_line.strip.split('\t').size.should eq(13) # "J" label + 12 PAF columns, no tags
      j_line.should_not match(/tp:A:P/)
    end

    it "prints full untruncated PAF rows (incl. trailing tags) with -p (bug fix)" do
      # Regression test for a real bug found in this port: the -p listing
      # must print the FULL row (JS: a[i].join("\t")), including trailing
      # SAM/PAF tags like tp:A:P and cm:i:, not truncate to the first 12 PAF
      # columns the way the -e error-report lines correctly do.
      code, out, _ = run_paftools("misjoin", "-p", "#{FIXTURES}/misjoin.paf")
      code.should eq(0)
      q1_line = out.lines.find { |l| l.starts_with?("q1\t5000000\t0\t1500000") }.not_nil!
      q1_line.strip.split('\t').size.should eq(14) # 12 PAF columns + 2 trailing tags
      q1_line.should match(/tp:A:P\tcm:i:100$/)
    end
  end

  # ── mapeval ──────────────────────────────────────────────────────────────
  describe "mapeval" do
    # 4 simulated reads (names encode ground truth), 3 distinct mapQ values

    it "buckets accuracy by mapping quality" do
      code, out, _ = run_paftools("mapeval", "#{FIXTURES}/mapeval.paf")
      code.should eq(0)
      out.should match(/^Q\t60\t3\t0\t0\.000000000\t3$/m)
      out.should match(/^Q\t0\t1\t1\t0\.250000000\t4$/m)
    end
  end

  # ── pafcmp ───────────────────────────────────────────────────────────────
  describe "pafcmp" do
    # base has r1(match),r2(will move),r3(missing from test); test has
    # r1(match),r2(moved→wrong),r4(new, not in base)

    it "reports wrong/missing base alignments" do
      code, out, _ = run_paftools("pafcmp", "#{FIXTURES}/pafcmp_base.paf", "#{FIXTURES}/pafcmp_test.paf")
      code.should eq(0)
      out.should match(/^W\tr2\t/m)
      out.should match(/^M\tr3\t/m)
      out.should match(/1 wrong test alignment/)
      out.should match(/1 base alignments missing/)
    end

    it "counts additional test alignments correctly (bug fix)" do
      # JS increments opt.n_out_high/n_out_low, but those fields only exist
      # on `eval`, not `opt` — a typo that makes the printed count always 0.
      # We track real counters instead, so r4 (new, mapQ>=10) must count as 1.
      _, out, _ = run_paftools("pafcmp", "#{FIXTURES}/pafcmp_base.paf", "#{FIXTURES}/pafcmp_test.paf")
      out.should match(/1 additional test alignments with mapQ>=10/)
    end
  end

  # ── ov-eval ──────────────────────────────────────────────────────────────
  describe "ov-eval" do
    # 3 to-ref alignments (r1,r3,r4 mutually overlapping); overlapper only
    # reports r1-r3, so r1-r4 and r3-r4 are inferred-but-missed

    it "reports overlap sensitivity" do
      code, out, _ = run_paftools("ov-eval", "-l", "400", "#{FIXTURES}/ov_toref.paf", "#{FIXTURES}/ov_ovlp.paf")
      code.should eq(0)
      out.should match(/^2 overlaps inferred from the reference mapping$/m)
      out.should match(/^1 missed by the read overlapper$/m)
      out.should match(/^50\.00% sensitivity$/m)
    end
  end

  # ── junceval / exoneval ──────────────────────────────────────────────────
  describe "junceval" do
    # 3-exon transcript T1 on chr1+; SAM read spans all 3 exons via 2 introns
    # that exactly match the annotated junctions

    it "reports 100% correct introns for an exact match" do
      code, out, _ = run_paftools("junceval", "#{FIXTURES}/junceval.gtf", "#{FIXTURES}/junceval.sam")
      code.should eq(0)
      out.should match(/# predicted introns: 2/)
      out.should match(/# correct introns: 2 \(100\.00%\)/)
    end
  end

  describe "exoneval" do
    it "reports 100% correct exons for an exact match" do
      code, out, _ = run_paftools("exoneval", "#{FIXTURES}/junceval.gtf", "#{FIXTURES}/junceval.sam")
      code.should eq(0)
      out.should match(/# predicted exons: 3/)
      out.should match(/# correct exons: 3 \(100\.00%\)/)
    end
  end

  # ── vcfstat ──────────────────────────────────────────────────────────────
  describe "vcfstat" do
    it "tabulates substitutions and indels" do
      code, out, _ = run_paftools("vcfstat", "#{FIXTURES}/vcfstat1.vcf")
      code.should eq(0)
      out.should match(/# substitutions: 3/)
      out.should match(/ts\/tv: 2\.000/)
      out.should match(/# insertions: 1/)
      out.should match(/# deletions: 1/)
    end

    it "treats a literal ALT ending in '>' as a real substitution, not a symbolic allele (bug fix)" do
      # JS's check is `a[0]=='<' || a[1]=='>'` — almost certainly meant to
      # test whether the ALT starts with '<' AND ends with '>' (a genuine
      # symbolic allele like <DEL>). As written, any two-char ALT whose
      # SECOND character happens to be '>' is also (wrongly) skipped. We use
      # a regex anchored on both ends instead, so "T>" is counted normally.
      code, out, _ = run_paftools("vcfstat", "#{FIXTURES}/vcfstat2.vcf")
      code.should eq(0)
      out.should match(/# substitutions: 2/)
    end
  end

  # ── vcfsel ───────────────────────────────────────────────────────────────
  describe "vcfsel" do
    it "filters VCF records by allele length" do
      code, out, _ = run_paftools("vcfsel", "-L", "20", "#{FIXTURES}/vcfsel.vcf")
      code.should eq(0)
      lines = out.strip.split('\n')
      lines.size.should eq(3) # header + 2 short records; the 50bp INS is excluded
    end
  end

  # ── sveval ───────────────────────────────────────────────────────────────
  describe "sveval" do
    it "scores SV calls against a truth VCF" do
      code, out, _ = run_paftools("sveval", "#{FIXTURES}/sveval_base.vcf", "#{FIXTURES}/sveval_call.vcf")
      code.should eq(0)
      out.should match(/^SN\t2\t1\t0\.500000$/m)
      out.should match(/^PC\t2\t1\t0\.500000$/m)
      out.should match(/^F1\t0\.500000$/m)
    end
  end

  # ── paf2gff ──────────────────────────────────────────────────────────────
  describe "paf2gff" do
    # spliced PAF: 500M100N500M100N500M → 3-exon transcript, identity=0.8529

    it "converts a spliced PAF alignment to GFF3" do
      code, out, _ = run_paftools("paf2gff", "#{FIXTURES}/paf2gff.paf")
      code.should eq(0)
      lines = out.strip.split('\n')
      lines[0].should match(/\ttranscript\t1001\t2700\t/)
      lines[0].should match(/identity=0\.8529/)
      lines.count { |l| l.includes?("\tCDS\t") }.should eq(3)
    end
  end

  # ── mason2fq ─────────────────────────────────────────────────────────────
  describe "mason2fq" do
    # simulated.1 forward @ chr1:1000 (1-based) → FASTQ header chr1!999!1099!+
    # simulated.2 reverse (flag 16) @ chr1:2000 → header chr1!1999!2099!- with
    # revcomp'd sequence and reversed quality string

    it "converts a mason2 SAM to FASTQ with correct revcomp on reverse reads" do
      code, out, _ = run_paftools("mason2fq", "#{FIXTURES}/mason2fq.sam")
      code.should eq(0)
      lines = out.strip.split('\n')
      lines[0].should eq("@1!chr1!999!1099!+ 0:0:0")
      lines[4].should eq("@2!chr1!1999!2099!- 0:0:0")
      lines[5].should eq("GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")
    end
  end

  # ── sim2bed ──────────────────────────────────────────────────────────────
  describe "sim2bed" do
    it "converts mason paired-end and long-read simulated names to BED" do
      code, out, _ = run_paftools("sim2bed", "#{FIXTURES}/sim2bed.txt")
      code.should eq(0)
      lines = out.strip.split('\n')
      lines[0].should eq("chr1\t1000\t1100\t@read1!chr1!1000_2000!1100_2100!F1/1\t0\tF")
      lines[1].should eq("chr2\t5000\t5500\t@read2!chr2!5000!5500!+\t0\t+")
    end
  end

  # ── pbsim2fq ─────────────────────────────────────────────────────────────
  describe "pbsim2fq" do
    it "converts a pbsim MAF alignment to FASTA with correct revcomp" do
      code, out, _ = run_paftools("pbsim2fq", "#{FIXTURES}/pbsim_ref.fa.fai", "#{FIXTURES}/pbsim_rc.maf")
      code.should eq(0)
      lines = out.strip.split('\n')
      lines[0].should eq(">S1_1!chr1!1000!1008!-")
      lines[1].should eq("CGGGTTTT") # revcomp of AAAACCCG
    end
  end

  # ── badread2fa ───────────────────────────────────────────────────────────
  describe "badread2fa" do
    # read1 (+strand, forward), read2 (-strand → coord-flipped), read3
    # (chimera → discarded)

    it "converts Badread FASTQ headers to FASTA, discarding chimeras" do
      code, out, err = run_paftools("badread2fa", "#{FIXTURES}/badread_ref.fai", "#{FIXTURES}/badread.fastq")
      code.should eq(0)
      lines = out.strip.split('\n')
      lines[0].should eq(">S1!chr1!1000!1500!+\tri:f:98.5")
      lines[2].should eq(">S2!chr1!97500!98000!-\tri:f:97.2")
      err.should match(/WARNING: discarded 1 reads/)
    end
  end
end
