module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # bedcov — compute bases in BED regions covered by target BED
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  private def self.read_bed_ivs(fn : String, to_merge : Bool) : {Hash(String, Array({Int32, Int32})), Hash(String, Array(Int32))}
    ivs_h = Hash(String, Array({Int32, Int32})).new
    File.open(fn) do |file|
      file.each_line(chomp: true) do |line|
        t = line.split('\t'); next if t.size < 3
        chr = t[0]; bst = t[1].to_i; ben = t[2].to_i
        ivs_h[chr] ||= [] of {Int32, Int32}
        if t.size >= 12 && t[9] =~ /^\d+/
          n = t[9].to_i; sz = t[10].split(','); sts = t[11].split(',')
          n.times { |j| ivs_h[chr] << {bst + sts[j].to_i, bst + sts[j].to_i + sz[j].to_i} }
        else
          ivs_h[chr] << {bst, ben}
        end
      end
    end
    idx_h = Hash(String, Array(Int32)).new
    ivs_h.each do |chr, ivs|
      intv_sort(ivs); intv_merge(ivs) if to_merge
      idx_h[chr] = intv_build(ivs)
    end
    {ivs_h, idx_h}
  end

  def self.cmd_bedcov(args : Array(String)) : Int32
    print_len = false; to_merge = true; fn_excl : String? = nil
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-p"; print_len = true
      when "-d"; to_merge = false
      when "-e"; i += 1; fn_excl = args[i]
      when "-h", "--help"
        STDERR.puts "Usage: paftools bedcov [-pd] <regions.bed> <target.bed>"; return 0
      else rest << args[i]
      end
      i += 1
    end
    if rest.size < 2
      STDERR.puts "Usage: paftools bedcov [-pd] <regions.bed> <target.bed>"; return 1
    end

    # Optional exclusion BED (always merged, never deduped) — regions overlapping
    # this file are excluded from the target-file feature segments.
    excl_ivs, excl_idx = fn_excl ? read_bed_ivs(fn_excl, true) : {Hash(String, Array({Int32, Int32})).new, Hash(String, Array(Int32)).new}

    # Read regions BED and build interval index per chromosome
    reg_ivs, reg_idx = read_bed_ivs(rest[0], to_merge)

    tot_len = 0_i64; hit_len = 0_i64

    open_in(rest[1]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('\t'); next if t.size < 3
        chr = t[0]; bst = t[1].to_i; ben = t[2].to_i
        segs = [] of {Int32, Int32}
        if t.size >= 12 && t[9] =~ /^\d+/
          n = t[9].to_i; sz = t[10].split(','); sts = t[11].split(',')
          n.times { |j| segs << {bst + sts[j].to_i, bst + sts[j].to_i + sz[j].to_i} }
        else
          segs << {bst, ben}
        end

        if fn_excl && excl_ivs.has_key?(chr)
          eivs = excl_ivs[chr]; eidx = excl_idx[chr]
          segs = segs.reject { |seg| !intv_ovlp(eivs, eidx, seg[0], seg[1]).empty? }
        end

        feat_len = segs.sum { |seg| seg[1] - seg[0] }.to_i64
        tot_len += feat_len

        next unless reg_ivs.has_key?(chr)
        ivs = reg_ivs[chr]; idx = reg_idx[chr]
        overlap = [] of {Int32, Int32}
        segs.each do |seg|
          intv_ovlp(ivs, idx, seg[0], seg[1]).each do |ovlp|
            overlap << {[ovlp[0], seg[0]].max, [ovlp[1], seg[1]].min}
          end
        end
        feat_hit = overlap.empty? ? 0_i64 : cov_len(overlap)
        hit_len += feat_hit
        puts (["F", t[0..3].join('\t'), feat_len, feat_hit]).join('\t') if print_len
      end
    end

    STDERR.puts "# target bases: #{tot_len}"
    STDERR.puts "# target bases overlapping regions: #{hit_len} (#{"%.2f" % (tot_len > 0 ? 100.0*hit_len/tot_len : 0.0)}%)"
    0
  end
end
