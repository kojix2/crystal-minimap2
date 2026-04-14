module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # liftover — lift BED coordinates via PAF CIGAR
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_liftover(args : Array(String)) : Int32
    to_merge = false; min_mapq = 5; min_len = 50000; max_div = 2.0
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-m"; to_merge = true
      when "-q"; i += 1; min_mapq = args[i].to_i
      when "-l"; i += 1; min_len = args[i].to_i
      when "-d"; i += 1; max_div = args[i].to_f
      when "-h", "--help"
        STDERR.puts "Usage: paftools liftover [-m] [-q INT] [-l INT] [-d FLOAT] <aln.paf> <query.bed>"
        return 0
      else rest << args[i]
      end
      i += 1
    end
    if rest.size < 2
      STDERR.puts "Usage: paftools liftover [-m] [-q INT] [-l INT] [-d FLOAT] <aln.paf> <query.bed>"
      return 1
    end

    # Read BED into per-chrom sorted (and optionally merged) interval lists
    bed = Hash(String, Array({Int32, Int32})).new
    bed_idx = Hash(String, Array(Int32)).new
    open_in(rest[1]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('\t'); next if t.size < 3
        chr = t[0]
        bed[chr] ||= [] of {Int32, Int32}
        bed[chr] << {t[1].to_i, t[2].to_i}
      end
    end
    bed.each do |chr, ivs|
      intv_sort(ivs)
      intv_merge(ivs) if to_merge
      bed_idx[chr] = intv_build(ivs)
    end

    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('\t')
        next if t.size < 12
        ivs = bed[t[0]]?; next unless ivs
        tp : String? = nil; cg : String? = nil
        t[12..].each do |tag|
          if m = /^(\S\S):[AZif]:(\S+)/.match(tag)
            tp = m[2] if m[1] == "tp"
            cg = m[2] if m[1] == "cg"
          end
        end
        next unless tp == "P" || tp == "I"
        next unless cg
        (1..3).each { |j| t[j] = t[j].to_i.to_s }
        (6..11).each { |j| t[j] = t[j].to_i.to_s }
        next if t[11].to_i < min_mapq || t[10].to_i < min_len

        regs = intv_ovlp(ivs, bed_idx[t[0]], t[2].to_i, t[3].to_i)
        next if regs.empty?

        if max_div >= 0.0 && max_div < 1.0
          n_gaps = 0; n_opens = 0
          cg.scan(/(\d+)([MID])/) do |mat|
            next if mat[2] == "M"
            n_gaps += mat[1].to_i; n_opens += 1
          end
          n_mm = t[10].to_i - t[9].to_i - n_gaps
          n_diff2 = n_mm + n_opens
          next if n_diff2.to_f / (n_diff2 + t[9].to_i) > max_div
        end

        strand = t[4]
        # Build sorted event list: each region contributes start and end events
        a = [] of {Int32, Int32, Int32, Int32} # {query_pos, type(0=st,1=en), reg_idx, ref_pos(-2=unset)}
        r = Array({Int32, Int32}).new(regs.size, {-2, -2})
        regs.each_with_index do |reg, reg_i|
          s = reg[0]; e = reg[1]
          if strand == "+"
            a << {s, 0, reg_i, -2}; a << {e - 1, 1, reg_i, -2}
          else
            a << {t[1].to_i - e, 0, reg_i, -2}; a << {t[1].to_i - s - 1, 1, reg_i, -2}
          end
        end
        a.sort_by! { |evnt| evnt[0] }

        # Walk CIGAR to map query positions to reference
        k = 0; rx = t[7].to_i
        y = strand == "+" ? t[2].to_i : t[1].to_i - t[3].to_i
        updated_a = a.map { |evnt| [evnt[0], evnt[1], evnt[2], evnt[3]] }

        cg.scan(/(\d+)([MID])/) do |mat|
          len = mat[1].to_i; op = mat[2][0]
          if op == 'D'
            rx += len; next
          end
          while k < updated_a.size && updated_a[k][0] < y
            k += 1
          end
          j = k
          while j < updated_a.size
            break if updated_a[j][0] >= y + len
            if y <= updated_a[j][0] && updated_a[j][0] < y + len
              updated_a[j][3] = op == 'M' ? rx + (updated_a[j][0] - y) : rx
            end
            j += 1
          end
          y += len; rx += len if op == 'M'
        end

        updated_a.each do |evnt|
          ri = evnt[2]
          if evnt[1] == 0
            r[ri] = {evnt[3], r[ri][1]}
          else
            r[ri] = {r[ri][0], evnt[3] + 1}
          end
        end
        r.each_with_index do |coords, reg_i|
          name = "#{t[0]}_#{regs[reg_i][0]}_#{regs[reg_i][1]}"
          r0 = coords[0]; r1 = coords[1]
          name += "_t5" if r0 < 0; r0 = t[7].to_i if r0 < 0
          name += "_t3" if r1 < 0; r1 = t[8].to_i if r1 < 0
          puts [t[5], r0, r1, name, 0, strand].join('\t')
        end
      end
    end
    0
  end
end
