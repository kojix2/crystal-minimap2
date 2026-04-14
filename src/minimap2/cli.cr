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
    out_file : String? = nil
    index_file : String? = nil # -d
    out_sam = false            # -a
    with_cigar = false         # -c
    with_eqx = false           # --eqx
    n_threads = 3
    show_version = false
    print_help = false

    parser = OptionParser.new do |p|
      p.summary_width = 12
      p.banner = "Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]"

      # ── Preset ──────────────────────────────────────────────────────────────
      p.on("-x STR", "Preset (map-ont, map-pb, map-hifi, asm5, asm10, asm20, " \
                     "splice, sr, ava-ont, ava-pb, ...)") do |v|
        preset = v
      end

      # ── Indexing ─────────────────────────────────────────────────────────────
      p.on("-H", "Use homopolymer-compressed (HPC) k-mers") do
        io.flag |= I_HPC
      end
      p.on("-k INT", "K-mer size [#{io.k}]") do |v|
        io.k = v.to_i
      end
      p.on("-w INT", "Minimizer window size [#{io.w}]") do |v|
        io.w = v.to_i
      end
      p.on("-I NUM", "Split index for every ~NUM input bases (k/m/g suffix) [8G]") do |v|
        io.batch_size = parse_num(v)
      end
      p.on("-d FILE", "Dump index to FILE") do |v|
        index_file = v
      end

      # ── Mapping ──────────────────────────────────────────────────────────────
      p.on("-f FLOAT", "Filter top FLOAT fraction of repetitive minimizers [#{mo.mid_occ_frac}]") do |v|
        mo.mid_occ_frac = v.to_f32
      end
      p.on("-g NUM", "Stop chain elongation if gap exceeds NUM bp [#{mo.max_gap}]") do |v|
        mo.max_gap = parse_num(v).to_i32
      end
      p.on("-r NUM", "Chaining bandwidth [#{mo.bw}]") do |v|
        parts = v.split(',')
        mo.bw = parts[0].to_i
        mo.bw_long = parts[1]?.try(&.to_i) || mo.bw_long
      end
      p.on("-n INT", "Minimal number of minimizers on a chain [#{mo.min_cnt}]") do |v|
        mo.min_cnt = v.to_i
      end
      p.on("-m INT", "Minimal chaining score [#{mo.min_chain_score}]") do |v|
        mo.min_chain_score = v.to_i
      end
      p.on("-X", "Skip self and dual mappings (all-vs-all mode)") do
        mo.flag |= F_NO_DUAL | F_NO_DIAG
      end
      p.on("-p FLOAT", "Min secondary-to-primary score ratio [#{mo.pri_ratio}]") do |v|
        mo.pri_ratio = v.to_f32
      end
      p.on("-N INT", "Retain at most INT secondary alignments [#{mo.best_n}]") do |v|
        mo.best_n = v.to_i
      end

      # ── Alignment ────────────────────────────────────────────────────────────
      p.on("-A INT", "Match score [#{mo.a}]") do |v|
        mo.a = v.to_i
      end
      p.on("-B INT", "Mismatch penalty [#{mo.b}]") do |v|
        mo.b = v.to_i
      end
      p.on("-O INT", "Gap open penalty; comma-sep for dual-affine [#{mo.q},#{mo.q2}]") do |v|
        parts = v.split(',')
        mo.q = parts[0].to_i
        mo.q2 = parts[1]?.try(&.to_i) || mo.q
      end
      p.on("-E INT", "Gap extend penalty; comma-sep for dual-affine [#{mo.e},#{mo.e2}]") do |v|
        parts = v.split(',')
        mo.e = parts[0].to_i
        mo.e2 = parts[1]?.try(&.to_i) || mo.e
      end
      p.on("-z INT", "Z-drop and inversion Z-drop [#{mo.zdrop},#{mo.zdrop_inv}]") do |v|
        parts = v.split(',')
        mo.zdrop = parts[0].to_i
        mo.zdrop_inv = parts[1]?.try(&.to_i) || mo.zdrop
      end
      p.on("-s INT", "Minimal peak DP alignment score [#{mo.min_dp_max}]") do |v|
        mo.min_dp_max = v.to_i
      end

      # ── Input / Output ───────────────────────────────────────────────────────
      p.on("-a", "Output in SAM format (PAF by default)") do
        out_sam = true
        mo.flag |= F_CIGAR
      end
      p.on("-o FILE", "Output alignments to FILE [stdout]") do |v|
        out_file = v
      end
      p.on("-c", "Output CIGAR in PAF") do
        with_cigar = true
        mo.flag |= F_CIGAR
      end
      p.on("--cs", "Output cs tag (short form)") do
        mo.flag |= F_OUT_CS
      end
      p.on("--eqx", "Write =/X CIGAR operators") do
        with_eqx = true
        mo.flag |= F_EQX
      end
      p.on("-t INT", "Number of threads [#{n_threads}]") do |v|
        n_threads = v.to_i
      end
      p.on("-K NUM", "Mini-batch size for mapping [500M]") do |v|
        mo.mini_batch_size = parse_num(v).to_i64
      end

      # ── Meta ─────────────────────────────────────────────────────────────────
      p.on("--version", "Show version and exit") do
        show_version = true
      end
      p.on("-h", "--help", "Show this help") do
        print_help = true
      end

      p.invalid_option do |flag|
        STDERR.puts "Error: unknown option #{flag}"
        STDERR.puts p
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

    # ── Open output ──────────────────────────────────────────────────────────
    out_io : IO = if (fn = out_file)
      File.open(fn, "w")
    else
      STDOUT
    end

    # ── Build / load index ───────────────────────────────────────────────────
    reader = MmIdxReader.new(target_fn, io)
    indices = [] of MmIdx
    while (mi = reader.read(n_threads))
      mapopt_update(mo, mi)
      if (idx_fn = index_file)
        File.open(idx_fn, "wb") { |f| mi.dump(f) }
      end
      indices << mi
    end
    reader.close

    if indices.empty?
      STDERR.puts "[ERROR] failed to load/build index from '#{target_fn}'"
      out_io.close unless out_file.nil?
      return 1
    end

    # Index-only mode (no queries).
    if query_fns.empty?
      out_io.close unless out_file.nil?
      return 0
    end

    # ── Write SAM header if -a ────────────────────────────────────────────────
    if out_sam
      write_sam_hdr(out_io, indices.first)
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
        seqs = bf.read_seqs(mo.mini_batch_size,
          with_qual: out_sam,
          with_comment: false,
          frag_mode: false)
        break if seqs.empty?

        # Pre-allocate result slots (one per read × per index).
        # We map every (seq, index) pair independently then write in order.
        n_pairs = seqs.size * indices.size
        # results[i] holds the alignments for pair i; nil means not yet done.
        results = Array(Array(MmReg1)?).new(n_pairs, nil)
        pending = Atomic(Int32).new(n_pairs)
        done_ch = Channel(Nil).new(1)

        seqs.each_with_index do |t, si|
          indices.each_with_index do |mi, ii|
            pair_idx = si * indices.size + ii
            map_ctx.spawn do
              results[pair_idx] = Minimap2.map(mi, t.l_seq, t.seq, mo, t.name)
              done_ch.send(nil) if pending.sub(1, :sequentially_consistent) == 1
            end
          end
        end

        done_ch.receive

        # Write results in original input order (deterministic output).
        seqs.each_with_index do |t, si|
          indices.each_with_index do |mi, ii|
            regs = results[si * indices.size + ii] || [] of MmReg1
            regs.each do |r|
              if out_sam
                write_sam(out_io, mi, t, r, regs.size, regs, mo.flag)
              else
                write_paf(out_io, mi, t, r, mo.flag)
              end
            end
            write_sam_unmapped(out_io, t) if out_sam && regs.empty?
          end
        end
      end

      bf.close
    end

    out_io.close unless out_file.nil?
    0
  end

  # ---------------------------------------------------------------------------
  # Write an unmapped SAM record (flag 4).
  # ---------------------------------------------------------------------------
  private def self.write_sam_unmapped(io : IO, t : BSeq1) : Nil
    io.print t.name
    io.print "\t4\t*\t0\t0\t*\t*\t0\t0\t"
    io.print t.seq
    io.print "\t"
    io.print t.qual || "*"
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
