require "option_parser"

module Minimap2
  # ---------------------------------------------------------------------------
  # Command-line interface — mirrors the minimap2 main() options
  # ---------------------------------------------------------------------------
  def self.run_cli(argv : Array(String)) : Int32
    # ── option state ──────────────────────────────────────────────────────────
    io = MmIdxOpt.new
    mo = MmMapOpt.new
    idxopt_init(io)
    mapopt_init(mo)

    preset : String? = nil
    argv.each_with_index do |arg, i|
      if arg == "-x"
        preset = argv[i + 1]?
      elsif arg.starts_with?("-x") && arg.size > 2
        preset = arg[2..]
      end
    end
    out_file : String? = nil
    index_file : String? = nil # -d
    rg_line : String? = nil    # -R
    out_sam = false            # -a
    with_cigar = false         # -c
    with_eqx = false           # --eqx
    n_threads = 3
    show_version = false
    print_help = false

    parser = OptionParser.new do |opt|
      opt.summary_width = 12
      opt.banner = "Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]"

      # ── Preset ──────────────────────────────────────────────────────────────
      opt.on("-x STR", "Preset (map-ont, map-pb, map-hifi, asm5, asm10, asm20, " \
                       "splice, sr, ava-ont, ava-pb, ...)") do |v|
        preset = v
      end

      # ── Indexing ─────────────────────────────────────────────────────────────
      opt.on("-H", "Use homopolymer-compressed (HPC) k-mers") do
        io.flag |= I_HPC
      end
      opt.on("-k INT", "K-mer size [#{io.k}]") do |v|
        io.k = v.to_i
      end
      opt.on("-w INT", "Minimizer window size [#{io.w}]") do |v|
        io.w = v.to_i
      end
      opt.on("-I NUM", "Split index for every ~NUM input bases (k/m/g suffix) [8G]") do |v|
        io.batch_size = parse_num(v)
      end
      opt.on("-d FILE", "Dump index to FILE") do |v|
        index_file = v
      end

      # ── Mapping ──────────────────────────────────────────────────────────────
      opt.on("-f FLOAT[,INT]", "Filter top FLOAT fraction of repetitive minimizers [#{mo.mid_occ_frac}]") do |v|
        # Mirrors main.c:325-331 — the optional second field sets max_occ.
        parts = v.split(',')
        value = parts[0].to_f64
        if value < 1.0
          mo.mid_occ_frac = value.to_f32
          mo.mid_occ = 0
        else
          mo.mid_occ = (value + 0.499).to_i32
        end
        if max_occ_str = parts[1]?
          mo.max_occ = max_occ_str.to_i
        end
      end
      opt.on("-g NUM", "Stop chain elongation if gap exceeds NUM bp [#{mo.max_gap}]") do |v|
        mo.max_gap = parse_num(v).to_i32
      end
      opt.on("-r NUM", "Chaining bandwidth [#{mo.bw}]") do |v|
        parts = v.split(',')
        mo.bw = parts[0].to_i
        mo.bw_long = parts[1]?.try(&.to_i) || mo.bw_long
      end
      opt.on("-n INT", "Minimal number of minimizers on a chain [#{mo.min_cnt}]") do |v|
        mo.min_cnt = v.to_i
      end
      opt.on("-m INT", "Minimal chaining score [#{mo.min_chain_score}]") do |v|
        mo.min_chain_score = v.to_i
      end
      opt.on("-X", "Skip self and dual mappings (all-vs-all mode)") do
        # Mirrors main.c: -X implies F_ALL_CHAINS | F_NO_DIAG | F_NO_DUAL | F_NO_LJOIN
        mo.flag |= F_ALL_CHAINS | F_NO_DIAG | F_NO_DUAL | F_NO_LJOIN
      end
      opt.on("-p FLOAT", "Min secondary-to-primary score ratio [#{mo.pri_ratio}]") do |v|
        mo.pri_ratio = v.to_f32
      end
      opt.on("-N INT", "Retain at most INT secondary alignments [#{mo.best_n}]") do |v|
        mo.best_n = v.to_i
      end

      # ── Alignment ────────────────────────────────────────────────────────────
      opt.on("-A INT", "Match score [#{mo.a}]") do |v|
        mo.a = v.to_i
      end
      opt.on("-B INT", "Mismatch penalty [#{mo.b}]") do |v|
        mo.b = v.to_i
      end
      opt.on("-O INT", "Gap open penalty; comma-sep for dual-affine [#{mo.q},#{mo.q2}]") do |v|
        parts = v.split(',')
        mo.q = parts[0].to_i
        mo.q2 = parts[1]?.try(&.to_i) || mo.q
      end
      opt.on("-E INT", "Gap extend penalty; comma-sep for dual-affine [#{mo.e},#{mo.e2}]") do |v|
        parts = v.split(',')
        mo.e = parts[0].to_i
        mo.e2 = parts[1]?.try(&.to_i) || mo.e
      end
      opt.on("-z INT", "Z-drop and inversion Z-drop [#{mo.zdrop},#{mo.zdrop_inv}]") do |v|
        parts = v.split(',')
        mo.zdrop = parts[0].to_i
        mo.zdrop_inv = parts[1]?.try(&.to_i) || mo.zdrop
      end
      opt.on("-s INT", "Minimal peak DP alignment score [#{mo.min_dp_max}]") do |v|
        mo.min_dp_max = v.to_i
      end
      opt.on("-u CHAR", "GT-AG direction: f=forward, b=both, n=none, r=reverse [b]") do |v|
        # Mirrors main.c:332-340 — splice presets start with both bits set,
        # so strand selection must clear the unwanted bit.
        case v
        when "f"
          mo.flag &= ~F_SPLICE_REV
        when "b"
          mo.flag |= F_SPLICE_FOR | F_SPLICE_REV
        when "n"
          mo.flag &= ~(F_SPLICE_FOR | F_SPLICE_REV)
        when "r"
          mo.flag &= ~F_SPLICE_FOR
        else
          STDERR.puts "[ERROR] unknown argument to -u: '#{v}' (expected f, b, n or r)"
          exit 1
        end
      end

      # ── Input / Output ───────────────────────────────────────────────────────
      opt.on("-a", "Output in SAM format (PAF by default)") do
        out_sam = true
        # Mirrors main.c: -a sets F_OUT_SAM | F_CIGAR
        mo.flag |= F_OUT_SAM | F_CIGAR
      end
      opt.on("-o FILE", "Output alignments to FILE [stdout]") do |v|
        out_file = v
      end
      opt.on("-c", "Output CIGAR in PAF") do
        with_cigar = true
        # Mirrors main.c: -c sets F_OUT_CG | F_CIGAR
        mo.flag |= F_OUT_CG | F_CIGAR
      end
      opt.on("--cs", "Output cs tag (short form)") do
        # Mirrors main.c: --cs implies CIGAR computation
        mo.flag |= F_OUT_CS | F_CIGAR
      end
      opt.on("--cs-long", "Output cs tag (long form)") do
        mo.flag |= F_OUT_CS | F_OUT_CS_LONG | F_CIGAR
      end
      opt.on("--MD", "Output MD tag") do
        # Mirrors main.c: --MD implies CIGAR computation
        mo.flag |= F_OUT_MD | F_CIGAR
      end
      opt.on("--eqx", "Write =/X CIGAR operators") do
        with_eqx = true
        mo.flag |= F_EQX
      end
      opt.on("-Y", "Use soft clipping for supplementary alignments") do
        mo.flag |= F_SOFTCLIP
      end
      opt.on("-L", "Write long CIGARs using the CG tag") do
        mo.flag |= F_LONG_CIGAR
      end
      opt.on("-R STR", "Read group header line (e.g. '@RG\\tID:foo\\tSM:bar')") do |v|
        rg_line = v
      end
      opt.on("-y", "Copy FASTA/Q comments to output") do
        mo.flag |= F_COPY_COMMENT
      end
      opt.on("--end-bonus INT", "Score bonus when alignment reaches end of read") do |v|
        mo.end_bonus = v.to_i
      end
      opt.on("-C INT", "Cost for a non-canonical splice site [#{mo.noncan}]") do |v|
        mo.noncan = v.to_i
      end
      opt.on("-G NUM", "Max intron length [#{mo.max_gap_ref}]") do |v|
        mo.max_gap_ref = parse_num(v).to_i32
      end
      opt.on("-F NUM", "Max fragment length [#{mo.max_frag_len}]") do |v|
        mo.max_frag_len = parse_num(v).to_i32
      end
      opt.on("--no-long-join", "Disable long join") do
        mo.flag |= F_NO_LJOIN
      end
      opt.on("--for-only", "Only map to forward strand of reference") do
        mo.flag |= F_FOR_ONLY
      end
      opt.on("--rev-only", "Only map to reverse strand of reference") do
        mo.flag |= F_REV_ONLY
      end
      opt.on("--no-end-flt", "Disable end filtering") do
        mo.flag |= F_NO_END_FLT
      end
      opt.on("--qstrand", "Segment query on the query strand") do
        mo.flag |= F_QSTRAND | F_NO_INV
      end
      opt.on("--secondary=yes|no", "Whether to output secondary alignments") do |v|
        # Mirrors main.c: "yes" clears the preset's F_NO_PRINT_2ND bit.
        case v
        when "no"
          mo.flag |= F_NO_PRINT_2ND
        when "yes"
          mo.flag &= ~F_NO_PRINT_2ND
        else
          STDERR.puts "[ERROR] unknown argument to --secondary: '#{v}' (expected yes or no)"
          exit 1
        end
      end
      opt.on("-t INT", "Number of threads [#{n_threads}]") do |v|
        n_threads = v.to_i
      end
      opt.on("-K NUM", "Mini-batch size for mapping [500M]") do |v|
        mo.mini_batch_size = parse_num(v).to_i64
      end
      opt.on("--split-prefix STR", "Prefix for split index") do |v|
        mo.split_prefix = v
      end
      opt.on("--dbg-seed-occ", "Print seed occurrence diagnostics") do
        @@dbg_flag |= DBG_SEED_FREQ
      end

      # ── Meta ─────────────────────────────────────────────────────────────────
      opt.on("--version", "Show version and exit") do
        show_version = true
      end
      opt.on("-h", "--help", "Show this help") do
        print_help = true
      end

      opt.invalid_option do |flag|
        STDERR.puts "Error: unknown option #{flag}"
        STDERR.puts opt
        exit 1
      end
    end

    # Apply preset BEFORE parsing remaining options so user can override.
    # (mirrors minimap2 main.c: mm_set_opt applied first then per-option overrides)
    if preset
      ret = set_opt(preset, io, mo)
      if ret < 0
        STDERR.puts "[ERROR] unknown preset '#{preset}'"
        return 1
      end
    end

    remaining = argv.dup
    begin
      parser.parse(remaining)
    rescue ex : OptionParser::InvalidOption
      STDERR.puts ex.message
      STDERR.puts parser
      return 1
    end

    if show_version
      puts "#{LIB_VERSION} (Crystal port)"
      return 0
    end

    if print_help || remaining.empty?
      puts parser
      return print_help ? 0 : 1
    end
    n_threads = 1 if n_threads < 1

    # ── Positional arguments ──────────────────────────────────────────────────
    target_fn = remaining.shift
    query_fns = remaining # may be empty (just build index)

    # Apply splice-flag cross-correction that doesn't need an index.
    if (mo.flag & F_SPLICE_FOR) != 0 || (mo.flag & F_SPLICE_REV) != 0
      mo.flag |= F_SPLICE
    end
    if check_opt(io, mo) < 0
      STDERR.puts "[ERROR] incompatible options"
      return 1
    end

    # ── Reject unimplemented feature combinations ─────────────────────────────
    # The audit identified these as silently accepted but producing wrong
    # output. Fail loudly until they are implemented.
    if (mo.flag & F_FRAG_MODE) != 0
      STDERR.puts "[ERROR] paired/fragment mapping is not yet output-compatible (see audit C-02/C-03/C-04)"
      return 1
    end
    if (mo.flag & F_SPLICE) != 0
      STDERR.puts "[ERROR] splice alignment is not yet output-compatible (see audit H-03)"
      return 1
    end
    if (mo.flag & F_SOFTCLIP) != 0
      STDERR.puts "[ERROR] -Y/--softclip is not yet implemented (see audit H-08)"
      return 1
    end
    if (mo.flag & F_LONG_CIGAR) != 0
      STDERR.puts "[ERROR] -L/--long-cigar is not yet implemented (see audit H-08)"
      return 1
    end
    if mo.split_prefix
      STDERR.puts "[ERROR] --split-prefix is not yet implemented (see audit C-09)"
      return 1
    end

    # ── Open output ──────────────────────────────────────────────────────────
    out_io : IO = if fn = out_file
      File.open(fn, "w")
    else
      STDOUT
    end

    # ── Build / load index ───────────────────────────────────────────────────
    # MmIdxReader keeps one output stream open, so multipart index dumps append
    # every part instead of repeatedly truncating the destination.
    reader = MmIdxReader.new(target_fn, io, index_file)
    indices = [] of MmIdx
    while mi = reader.read(n_threads)
      mapopt_update(mo, mi)
      indices << mi
    end
    reader.close

    if indices.empty?
      STDERR.puts "[ERROR] failed to load/build index from '#{target_fn}'"
      out_io.close unless out_file.nil?
      return 1
    end

    # Index-only mode: -d without query files is fine; otherwise missing input
    # is an error (mirrors main.c:430-433).
    if query_fns.empty?
      out_io.close unless out_file.nil?
      if index_file.nil? && MmIdx.idx?(target_fn) == 0_i64
        STDERR.puts "[ERROR] missing input: please specify a query file to map or use -d to dump the index"
        return 1
      end
      return 0
    end

    if indices.size > 1
      STDERR.puts "[ERROR] mapping against a multipart index is not yet supported (see audit C-09)"
      out_io.close unless out_file.nil?
      return 1
    end

    # ── Write SAM header if -a ────────────────────────────────────────────────
    rg_id : String? = nil
    if out_sam
      write_sam_hdr(out_io, indices.first, LIB_VERSION, rg_line)
      rg_id = Minimap2.rg_id(rg_line)
    end

    # ── Parallel mapping context ──────────────────────────────────────────────
    # One Parallel context shared across all query files; resized to n_threads.
    map_ctx = Fiber::ExecutionContext::Parallel.new("mapper", n_threads)

    # ── Map each query file ───────────────────────────────────────────────────
    query_fns.each do |qry_fn|
      bf = BSeqFile.open(qry_fn)
      unless bf
        STDERR.puts "[ERROR] cannot open query file '#{qry_fn}'"
        next
      end

      loop do
        is_frag = (mo.flag & F_FRAG_MODE) != 0
        seqs = bf.read_seqs(mo.mini_batch_size,
          with_qual: out_sam,
          with_comment: (mo.flag & F_COPY_COMMENT) != 0,
          frag_mode: is_frag)
        break if seqs.empty?

        if is_frag && query_fns.size == 1
          # Fragment/paired-end mode: group consecutive reads by name
          # and map as multi-segment fragments
          frag_groups = group_fragments(seqs)

          frag_groups.each_with_index do |group, grp_i|
            indices.each_with_index do |midx, mi_i|
              n_segs = group.size
              if n_segs == 1
                # Single-end: map normally
                bseq = group[0]
                tbuf = MmTbuf.new
                regs = Minimap2.map(midx, bseq.l_seq, bseq.seq, mo, bseq.name, tbuf)
                regs.each do |reg|
                  write_output(out_io, out_sam, midx, bseq, reg, regs.size, regs, mo, tbuf.rep_len, rg_id)
                end
                write_sam_unmapped(out_io, bseq, mo.flag, rg_id) if out_sam && regs.empty?
              else
                # Paired-end: apply pe_ori revcomp and map as fragment
                seqs_str = Array(String).new(n_segs)
                qlens = Array(Int32).new(n_segs)
                n_segs.times do |seg_i|
                  s = group[seg_i].seq
                  # pe_ori: bit 1 = revcomp read1, bit 0 = revcomp read2
                  if seg_i == 0 && (mo.pe_ori & 2) != 0
                    s = s.reverse.tr("ACGTacgtNn", "TGCAtgcaNn")
                  elsif seg_i == 1 && (mo.pe_ori & 1) != 0
                    s = s.reverse.tr("ACGTacgtNn", "TGCAtgcaNn")
                  end
                  seqs_str << s
                  qlens << group[seg_i].l_seq
                end

                tbuf = MmTbuf.new
                all_regs = Minimap2.map_frag(midx, seqs_str, qlens, mo, group[0].name, tbuf)

                n_segs.times do |seg_i|
                  bseq = group[seg_i]
                  regs = all_regs[seg_i]
                  regs.each do |reg|
                    write_output(out_io, out_sam, midx, bseq, reg, regs.size, regs, mo, tbuf.rep_len, rg_id)
                  end
                  write_sam_unmapped(out_io, bseq, mo.flag, rg_id) if out_sam && regs.empty?
                end
              end
            end
          end
        else
          # Independent mapping of each sequence
          results, rep_lens = map_batch(map_ctx, indices, seqs, mo, n_threads)

          seqs.each_with_index do |bseq, seq_i|
            indices.each_with_index do |midx, mi_i|
              pair_idx = seq_i * indices.size + mi_i
              regs = results[pair_idx] || [] of MmReg1
              regs.each do |reg|
                write_output(out_io, out_sam, midx, bseq, reg, regs.size, regs, mo, rep_lens[pair_idx], rg_id)
              end
              write_sam_unmapped(out_io, bseq, mo.flag, rg_id) if out_sam && regs.empty?
            end
          end
        end # else (not frag mode)
      end   # loop do

      bf.close
    end

    out_io.close unless out_file.nil?
    0
  end

  # Map a mini-batch while keeping output ordering deterministic.  The parallel
  # path uses a fixed worker pool instead of spawning one fiber per read/index
  # pair, which keeps scheduling overhead bounded for large batches.
  private def self.map_batch(map_ctx : Fiber::ExecutionContext::Parallel,
                             indices : Array(MmIdx), seqs : Array(BSeq1),
                             mo : MmMapOpt, n_threads : Int32) : {Array(Array(MmReg1)?), Array(Int32)}
    n_pairs = seqs.size * indices.size
    results = Array(Array(MmReg1)?).new(n_pairs, nil)
    rep_lens = Array(Int32).new(n_pairs, 0)
    return {results, rep_lens} if n_pairs == 0

    if n_threads <= 1 || n_pairs == 1
      n_pairs.times do |pair_idx|
        map_one_pair(indices, seqs, mo, pair_idx, results, rep_lens)
      end
      return {results, rep_lens}
    end

    worker_count = {n_threads, n_pairs}.min
    next_pair = Atomic(Int32).new(0)
    pending = Atomic(Int32).new(worker_count)
    done_ch = Channel(Nil).new(1)

    worker_count.times do
      map_ctx.spawn do
        loop do
          pair_idx = next_pair.add(1, :relaxed)
          break if pair_idx >= n_pairs
          map_one_pair(indices, seqs, mo, pair_idx, results, rep_lens)
        end
        done_ch.send(nil) if pending.sub(1, :acquire_release) == 1
      end
    end

    done_ch.receive
    {results, rep_lens}
  end

  private def self.map_one_pair(indices : Array(MmIdx), seqs : Array(BSeq1),
                                mo : MmMapOpt, pair_idx : Int32,
                                results : Array(Array(MmReg1)?),
                                rep_lens : Array(Int32)) : Nil
    mi_i = pair_idx % indices.size
    seq_i = pair_idx // indices.size
    midx = indices[mi_i]
    bseq = seqs[seq_i]

    # Do not swallow mapping errors: an internal failure (index corruption,
    # range error, DP bug) must not be silently converted into an "unmapped"
    # read with a successful exit status (see audit H-09).
    tbuf = MmTbuf.new
    results[pair_idx] = Minimap2.map(midx, bseq.l_seq, bseq.seq, mo, bseq.name, tbuf)
    rep_lens[pair_idx] = tbuf.rep_len
  end

  # ---------------------------------------------------------------------------
  # Write a single alignment result to the output.
  # ---------------------------------------------------------------------------
  private def self.write_output(out_io : IO, out_sam : Bool,
                                midx : MmIdx, bseq : BSeq1,
                                reg : MmReg1, n_regs : Int32,
                                all_regs : Array(MmReg1),
                                mo : MmMapOpt, rep_len : Int32,
                                rg_id : String? = nil) : Nil
                  return if (mo.flag & F_NO_PRINT_2ND) != 0 && reg.parent >= 0 && reg.parent != reg.id
    need_seq = (mo.flag & (F_OUT_CS | F_OUT_CS_LONG | F_OUT_MD)) != 0
    if out_sam
      if need_seq && reg.p
        tseq_arr = get_ref_seq(midx, reg.rid, reg.rs, reg.re, false)
        qseq_arr = encode_seq(bseq.seq[reg.qs...reg.qe])
        qseq_arr = rev_comp(qseq_arr) if reg.rev?
        write_sam(out_io, midx, bseq, reg, n_regs, all_regs, mo.flag, tseq_arr, qseq_arr, rg_id)
      else
        write_sam(out_io, midx, bseq, reg, n_regs, all_regs, mo.flag, nil, nil, rg_id)
      end
    else
      if need_seq && reg.p
        tseq_arr = get_ref_seq(midx, reg.rid, reg.rs, reg.re, false)
        qseq_arr = encode_seq(bseq.seq[reg.qs...reg.qe])
        qseq_arr = rev_comp(qseq_arr) if reg.rev?
        write_paf(out_io, midx, bseq, reg, mo.flag, rep_len, tseq_arr, qseq_arr)
      else
        write_paf(out_io, midx, bseq, reg, mo.flag, rep_len)
      end
    end
  end

  # ---------------------------------------------------------------------------
  # Group consecutive reads by name (paired-end / fragment mode).
  # Reads with the same base name (ignoring /1 /2 suffix) are grouped.
  # ---------------------------------------------------------------------------
  private def self.group_fragments(seqs : Array(BSeq1)) : Array(Array(BSeq1))
    groups = [] of Array(BSeq1)
    return groups if seqs.empty?

    current = [seqs[0]]
    (1...seqs.size).each do |i|
      if same_qname_base?(seqs[i].name, current[0].name)
        current << seqs[i]
      else
        groups << current
        current = [seqs[i]]
      end
    end
    groups << current
    groups
  end

  # Check if two names match ignoring /1 /2 suffix.
  private def self.same_qname_base?(a : String, b : String) : Bool
    return false if a.empty? || b.empty?
    la = a.size; lb = b.size
    ea = (la >= 2 && a[la - 2] == '/' && (a[la - 1] == '1' || a[la - 1] == '2')) ? la - 2 : la
    eb = (lb >= 2 && b[lb - 2] == '/' && (b[lb - 1] == '1' || b[lb - 1] == '2')) ? lb - 2 : lb
    return false if ea != eb
    a[0...ea] == b[0...eb]
  end

  # ---------------------------------------------------------------------------
  # Write an unmapped SAM record (flag 4).
  # ---------------------------------------------------------------------------
  private def self.write_sam_unmapped(io : IO, t : BSeq1, opt_flag : Int64,
                                      rg_id : String? = nil) : Nil
    io.print t.name
    io.print "\t4\t*\t0\t0\t*\t*\t0\t0\t"
    io.print t.seq
    io.print "\t"
    io.print (opt_flag & F_NO_QUAL) == 0 ? (t.qual || "*") : "*"
    if rg_id
      io.print "\tRG:Z:"; io.print rg_id
    end
    if (opt_flag & F_COPY_COMMENT) != 0 && (comment = t.comment)
      io.print "\t"; io.print comment
    end
    io.print "\n"
  end

  # ---------------------------------------------------------------------------
  # Parse a number with optional k/m/g suffix.
  # ---------------------------------------------------------------------------
  private def self.parse_num(s : String) : UInt64
    s = s.strip
    mult = case s[-1].downcase
           when 'k' then 1_000_u64
           when 'm' then 1_000_000_u64
           when 'g' then 1_000_000_000_u64
           else          1_u64
           end
    base = mult == 1_u64 ? s : s[0..-2]
    (base.to_f64 * mult).to_u64
  end
end
