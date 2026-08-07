require "option_parser"
require "compress/gzip"

# paftools.cr — Crystal port of paftools.js (minimap2 companion utilities)
module Paftools
  VERSION = "2.30 (Crystal port)"
end

require "./paftools/common"
require "./paftools/stat"
require "./paftools/view"
require "./paftools/convert"
require "./paftools/liftover"
require "./paftools/call"
require "./paftools/annotation"
require "./paftools/asm"
require "./paftools/bedcov"
require "./paftools/vcfpair"
require "./paftools/misjoin"
require "./paftools/mapeval"
require "./paftools/pafcmp"
require "./paftools/simreads"
require "./paftools/junceval"
require "./paftools/ov_eval"
require "./paftools/vcfextra"
require "./paftools/paf2gff"

module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # run — command dispatcher
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.run(argv : Array(String)) : Int32
    if argv.empty?
      STDERR.puts "Usage: paftools <command> [arguments]"
      STDERR.puts "Commands:"
      STDERR.puts "  view       convert PAF to BLAST-like (default), MAF, or LASTZ-cigar"
      STDERR.puts "  sam2paf    convert SAM to PAF"
      STDERR.puts "  delta2paf  convert MUMmer delta to PAF"
      STDERR.puts "  gff2bed    convert GFF/GTF to BED12"
      STDERR.puts "  gff2junc   convert GFF3 to junction BED"
      STDERR.puts "  splice2bed convert splice alignments (PAF/SAM) to BED12"
      STDERR.puts "  longcs2seq reconstruct sequences from long-cs PAF"
      STDERR.puts ""
      STDERR.puts "  stat       basic mapping statistics (PAF or SAM)"
      STDERR.puts "  asmstat    assembly alignment statistics"
      STDERR.puts "  asmgene    gene completeness evaluation"
      STDERR.puts "  liftover   lift BED coordinates via PAF CIGAR"
      STDERR.puts "  call       call variants from PAF with cs tag"
      STDERR.puts "  bedcov     compute coverage of BED regions"
      STDERR.puts "  vcfpair    pair diploid VCF haplotypes"
      STDERR.puts "  misjoin    detect candidate misjoins/inversions in an asm-to-ref PAF"
      STDERR.puts "  mapeval    evaluate mapping accuracy on simulated reads"
      STDERR.puts "  pafcmp     compare two PAF files (base vs. test)"
      STDERR.puts "  ov-eval    estimate read-overlap sensitivity from a to-ref mapping"
      STDERR.puts "  junceval   evaluate splice junction accuracy against a GTF"
      STDERR.puts "  exoneval   evaluate exon-level accuracy against a GTF"
      STDERR.puts "  vcfstat    tabulate substitution/indel counts from a VCF"
      STDERR.puts "  vcfsel     select VCF records by allele length"
      STDERR.puts "  sveval     evaluate structural-variant calls against truth VCF"
      STDERR.puts "  paf2gff    convert a spliced-alignment PAF to GFF3"
      STDERR.puts "  mason2fq   convert mason2-simulated SAM to FASTQ"
      STDERR.puts "  sim2bed    convert simulated-read names to BED"
      STDERR.puts "  pbsim2fq   convert pbsim MAF alignments to FASTA"
      STDERR.puts "  badread2fa convert Badread-simulated FASTQ/FASTA to FASTA"
      STDERR.puts ""
      STDERR.puts "  version    print version"
      return 1
    end
    cmd = argv.shift
    case cmd
    when "stat"                 then cmd_stat(argv)
    when "view"                 then cmd_view(argv)
    when "sam2paf"              then cmd_sam2paf(argv)
    when "delta2paf"            then cmd_delta2paf(argv)
    when "liftover", "liftOver" then cmd_liftover(argv)
    when "call"                 then cmd_call(argv)
    when "longcs2seq"           then cmd_longcs2seq(argv)
    when "gff2bed"              then cmd_gff2bed(argv)
    when "gff2junc"             then cmd_gff2junc(argv)
    when "splice2bed"           then cmd_splice2bed(argv)
    when "asmstat"              then cmd_asmstat(argv)
    when "asmgene"              then cmd_asmgene(argv)
    when "bedcov"               then cmd_bedcov(argv)
    when "vcfpair"              then cmd_vcfpair(argv)
    when "misjoin"              then cmd_misjoin(argv)
    when "mapeval"              then cmd_mapeval(argv)
    when "pafcmp"               then cmd_pafcmp(argv)
    when "mason2fq"             then cmd_mason2fq(argv)
    when "sim2bed"              then cmd_sim2bed(argv)
    when "pbsim2fq"             then cmd_pbsim2fq(argv)
    when "badread2fa"           then cmd_badread2fa(argv)
    when "junceval"             then cmd_junceval(argv)
    when "exoneval"             then cmd_exoneval(argv)
    when "ov-eval"              then cmd_ov_eval(argv)
    when "vcfstat"              then cmd_vcfstat(argv)
    when "vcfsel"               then cmd_vcfsel(argv)
    when "sveval"               then cmd_sveval(argv)
    when "paf2gff"              then cmd_paf2gff(argv)
    when "version"              then puts VERSION; 0
    else                             STDERR.puts "Error: unknown command '#{cmd}'"; 1
    end
  end
end
