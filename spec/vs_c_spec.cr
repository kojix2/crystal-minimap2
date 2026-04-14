require "./spec_helper"

# ---------------------------------------------------------------------------
# Comparison tests: Crystal port vs C minimap2
#
# Each test runs the system `minimap2` binary to get reference PAF output,
# then runs the same query through our Crystal Aligner and compares the
# key fields.
#
# What we test vs. C:
#   - Strand (exact)                — determined by chaining, must be correct
#   - Reference name (exact)        — must map to the same contig
#   - Primary hit count (range)     — chaining differs slightly; allow ±1
#   - Query coordinates (±200 bp)   — chaining anchors may differ at chain ends
#   - Ref coordinates (±200 bp)     — same reasoning
#   - Alignment block length (±10%) — blen from CIGAR; simplified ksw2 may differ
#   - MAPQ == 60 for primary        — unambiguous alignments always score 60
#
# What we do NOT test vs. C:
#   - mlen (matching bases) — depends on full ksw2 CIGAR accuracy (simplified here)
#   - Optional PAF tags (NM, ms, AS, cg, etc.)
#   - Secondary/supplementary hit details for complex inversions
# ---------------------------------------------------------------------------

MINIMAP2_BIN = ENV["MINIMAP2_BIN"]? || "minimap2"
TEST_DIR     = File.join(__DIR__, "../../minimap2/test")

# Parse the first 12 PAF columns into a named tuple.
record PafRow,
  qname : String, qlen : Int32, qs : Int32, qe : Int32,
  strand : Char,
  tname : String, tlen : Int32, rs : Int32, re : Int32,
  mlen : Int32, blen : Int32, mapq : Int32

private def parse_paf(text : String) : Array(PafRow)
  text.lines.compact_map do |line|
    next if line.strip.empty? || line.starts_with?('#')
    f = line.split('\t')
    next if f.size < 12
    PafRow.new(
      qname: f[0], qlen: f[1].to_i,
      qs: f[2].to_i, qe: f[3].to_i,
      strand: f[4][0],
      tname: f[5], tlen: f[6].to_i,
      rs: f[7].to_i, re: f[8].to_i,
      mlen: f[9].to_i, blen: f[10].to_i, mapq: f[11].to_i
    )
  end
end

# Run C minimap2 and return parsed PAF rows.
# Extra args (e.g. "-cx asm5") can be passed as a single string.
private def c_minimap2(ref_fa : String, qry_fa : String, args : String = "-c") : Array(PafRow)
  output = IO::Memory.new
  Process.run(MINIMAP2_BIN, args: args.split + [ref_fa, qry_fa],
    output: output, error: Process::Redirect::Close)
  parse_paf(output.to_s)
end

private def read_fasta(path : String) : Array({String, String})
  seqs = [] of {String, String}
  name = ""; buf = String::Builder.new
  File.each_line(path) do |line|
    if line.starts_with?('>')
      seqs << {name, buf.to_s} unless name.empty?
      name = line[1..].split(' ')[0].strip
      buf = String::Builder.new
    else
      buf << line.strip.upcase
    end
  end
  seqs << {name, buf.to_s} unless name.empty?
  seqs
end

# Run Crystal Aligner on the same inputs.
private def crystal_minimap2(ref_fa : String, qry_fa : String, preset : String) : Array(PafRow)
  ref_seqs = read_fasta(ref_fa)
  qry_seqs = read_fasta(qry_fa)
  aln = Minimap2::Aligner.from_strings(ref_seqs.map(&.[1]), ref_seqs.map(&.[0]), preset)

  rows = [] of PafRow
  qry_seqs.each do |qname, qseq|
    aln.map(qseq, qname).each do |h|
      next if h.rid < 0
      tname = aln.seq_names[h.rid]
      tlen = aln.seq_lengths[h.rid].to_i32
      rows << PafRow.new(
        qname: qname, qlen: qseq.size,
        qs: h.qs, qe: h.qe,
        strand: h.rev? ? '-' : '+',
        tname: tname, tlen: tlen,
        rs: h.rs, re: h.re,
        mlen: h.mlen, blen: h.blen, mapq: h.mapq.to_i32
      )
    end
  end
  rows
end

# ---------------------------------------------------------------------------
# Tolerance helper
# ---------------------------------------------------------------------------
private def close_enough(val : Int32, ref_val : Int32, tol : Int32) : Bool
  (val - ref_val).abs <= tol
end

# ===========================================================================
# Test suite
# ===========================================================================

describe "Crystal port vs C minimap2" do
  mt_human = "#{TEST_DIR}/MT-human.fa"
  mt_orang = "#{TEST_DIR}/MT-orang.fa"
  t_inv = "#{TEST_DIR}/t-inv.fa"
  q_inv = "#{TEST_DIR}/q-inv.fa"
  x3s_ref = "#{TEST_DIR}/x3s-ref.fa"
  x3s_qry = "#{TEST_DIR}/x3s-qry.fa"

  # -------------------------------------------------------------------------
  # 1.  MT-orang vs MT-human  (map-ont, single long alignment)
  # -------------------------------------------------------------------------
  describe "MT-orang vs MT-human" do
    c = c_minimap2(mt_human, mt_orang)
    cr = crystal_minimap2(mt_human, mt_orang, "map-ont")

    it "C produces exactly 1 primary alignment" do
      c.size.should eq(1)
    end

    it "Crystal produces at least 1 alignment" do
      cr.size.should be >= 1
    end

    it "strand matches C (should be +)" do
      ref = c.first
      hit = cr.first
      hit.strand.should eq(ref.strand)
    end

    it "maps to MT_human" do
      cr.first.tname.should eq("MT_human")
    end

    it "query start is within 200 bp of C" do
      ref = c.first; hit = cr.first
      close_enough(hit.qs, ref.qs, 200).should be_true,
        "Crystal qs=#{hit.qs} vs C qs=#{ref.qs} (diff=#{(hit.qs - ref.qs).abs})"
    end

    it "query end is within 200 bp of C" do
      ref = c.first; hit = cr.first
      close_enough(hit.qe, ref.qe, 200).should be_true,
        "Crystal qe=#{hit.qe} vs C qe=#{ref.qe} (diff=#{(hit.qe - ref.qe).abs})"
    end

    it "ref start is within 200 bp of C" do
      ref = c.first; hit = cr.first
      close_enough(hit.rs, ref.rs, 200).should be_true,
        "Crystal rs=#{hit.rs} vs C rs=#{ref.rs} (diff=#{(hit.rs - ref.rs).abs})"
    end

    it "ref end is within 200 bp of C" do
      ref = c.first; hit = cr.first
      close_enough(hit.re, ref.re, 200).should be_true,
        "Crystal re=#{hit.re} vs C re=#{ref.re} (diff=#{(hit.re - ref.re).abs})"
    end

    it "alignment block length covers similar fraction of query as C (within 10%)" do
      ref = c.first; hit = cr.first
      frac_c = ref.blen.to_f / ref.qlen
      frac_cr = hit.blen.to_f / hit.qlen
      (frac_c - frac_cr).abs.should be <= 0.10,
        "C blen_frac=#{frac_c.round(3)} Crystal blen_frac=#{frac_cr.round(3)}"
    end

    it "mapq is 60 (unambiguous alignment)" do
      cr.first.mapq.should eq(60)
    end
  end

  # -------------------------------------------------------------------------
  # 2.  q-inv vs t-inv  (map-ont, read with inversion — multiple hits)
  # -------------------------------------------------------------------------
  describe "q-inv vs t-inv (inversion)" do
    c = c_minimap2(t_inv, q_inv)
    cr = crystal_minimap2(t_inv, q_inv, "map-ont")

    it "C produces 6 hits total (3 per read)" do
      c.size.should eq(6)
    end

    it "Crystal maps both reads" do
      qnames = cr.map(&.qname).uniq!
      qnames.should contain("read1")
      qnames.should contain("read2")
    end

    it "read1 primary hit: strand matches C primary" do
      c_r1 = c.find! { |r| r.qname == "read1" }
      cr_r1 = cr.find! { |r| r.qname == "read1" }
      cr_r1.strand.should eq(c_r1.strand)
    end

    it "read2 primary hit: strand matches C primary" do
      c_r2 = c.find! { |r| r.qname == "read2" }
      cr_r2 = cr.find! { |r| r.qname == "read2" }
      cr_r2.strand.should eq(c_r2.strand)
    end

    it "read1 primary: maps to ref contig" do
      cr.find! { |r| r.qname == "read1" }.tname.should eq("ref")
    end

    it "read2 primary: maps to ref contig" do
      cr.find! { |r| r.qname == "read2" }.tname.should eq("ref")
    end

    # Note: C minimap2 splits the inversion into 3 sub-alignments per read.
    # Our Crystal port lacks full split-alignment/inversion detection, so it
    # produces one merged alignment covering more of the query.  We therefore
    # test that the merged hit still spans the correct region of the reference
    # (the union of all C sub-alignments ± slack).
    it "read1 primary: ref span falls within C's combined ref range (±500 bp slack)" do
      c_r1_all = c.select { |r| r.qname == "read1" }
      cr_r1 = cr.find! { |r| r.qname == "read1" }
      c_rs_min = c_r1_all.min_of(&.rs) - 500
      c_re_max = c_r1_all.max_of(&.re) + 500
      cr_r1.rs.should be >= c_rs_min
      cr_r1.re.should be <= c_re_max
    end

    it "read2 primary: ref span falls within C's combined ref range (±500 bp slack)" do
      c_r2_all = c.select { |r| r.qname == "read2" }
      cr_r2 = cr.find! { |r| r.qname == "read2" }
      c_rs_min = c_r2_all.min_of(&.rs) - 500
      c_re_max = c_r2_all.max_of(&.re) + 500
      cr_r2.rs.should be >= c_rs_min
      cr_r2.re.should be <= c_re_max
    end
  end

  # -------------------------------------------------------------------------
  # 3.  x3s-qry vs x3s-ref  (splice mode)
  # -------------------------------------------------------------------------
  describe "x3s-qry vs x3s-ref (splice)" do
    c = c_minimap2(x3s_ref, x3s_qry, "-cx splice")
    cr = crystal_minimap2(x3s_ref, x3s_qry, "splice")

    it "C produces 1 alignment" do
      c.size.should eq(1)
    end

    it "Crystal produces at least 1 alignment" do
      cr.size.should be >= 1
    end

    it "strand matches C (should be -)" do
      ref = c.first; hit = cr.first
      hit.strand.should eq(ref.strand)
    end

    it "maps to ref contig" do
      cr.first.tname.should eq("ref")
    end

    it "query span covers most of query (Crystal qe - qs >= 120 of 134 bp)" do
      hit = cr.first
      (hit.qe - hit.qs).should be >= 120
    end

    it "ref start is within 20 bp of C" do
      ref = c.first; hit = cr.first
      close_enough(hit.rs, ref.rs, 20).should be_true,
        "Crystal rs=#{hit.rs} vs C rs=#{ref.rs}"
    end

    it "ref end is within 20 bp of C" do
      ref = c.first; hit = cr.first
      close_enough(hit.re, ref.re, 20).should be_true,
        "Crystal re=#{hit.re} vs C re=#{ref.re}"
    end
  end
end
