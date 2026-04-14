module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # bedcov — compute bases in BED regions covered by target BED
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_bedcov(args : Array(String)) : Int32
    print_len = false; to_merge = true
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-p"; print_len = true
      when "-d"; to_merge = false
      when "-h", "--help"
        STDERR.puts "Usage: paftools bedcov [-pd] <regions.bed> <target.bed>"; return 0
      else rest << args[i]
      end
      i += 1
    end
    if rest.size < 2
      STDERR.puts "Usage: paftools bedcov [-pd] <regions.bed> <target.bed>"; return 1
    end

    # Read regions BED and build interval index per chromosome
    reg_ivs = Hash(String, Array({Int32, Int32})).new
    reg_idx = Hash(String, Array(Int32)).new
    File.open(rest[0]) do |file|
      file.each_line(chomp: true) do |line|
        t = line.split('\t'); next if t.size < 3
        chr = t[0]; bst = t[1].to_i; ben = t[2].to_i
        reg_ivs[chr] ||= [] of {Int32, Int32}
        if t.size >= 12 && t[9] =~ /^\d+/
          n = t[9].to_i; sz = t[10].split(','); sts = t[11].split(',')
          n.times { |j| reg_ivs[chr] << {bst + sts[j].to_i, bst + sts[j].to_i + sz[j].to_i} }
        else
          reg_ivs[chr] << {bst, ben}
        end
      end
    end
    reg_ivs.each do |chr, ivs|
      intv_sort(ivs); intv_merge(ivs) if to_merge
      reg_idx[chr] = intv_build(ivs)
    end

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
