module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # vcfpair — pair diploid VCF haplotypes
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_vcfpair(args : Array(String)) : Int32
    is_male = false; sample = "syndip"; hgver : String? = nil
    par = {"37" => [{0, 2699520}, {154931043, 155260560}]}
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-m"; is_male = true
      when "-s"; i += 1; sample = args[i]
      when "-g"; i += 1; hgver = args[i]
      when "-h", "--help"
        STDERR.puts "Usage: paftools vcfpair [-m] [-s sample] [-g 37] <in.pair.vcf>"; return 0
      else rest << args[i]
      end
      i += 1
    end
    if rest.empty?
      STDERR.puts "Usage: paftools vcfpair [-m] [-s sample] [-g 37] <in.pair.vcf>"; return 1
    end
    re_ctg = is_male ? /^(chr)?([0-9]+|X|Y)$/ : /^(chr)?([0-9]+|X)$/

    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        if line.starts_with?('#')
          next if line =~ /^##(source|reference)=/
          if m = /^##contig=.*ID=([^\s,]+)/.match(line)
            next unless re_ctg.match(m[1])
          elsif line.starts_with?("#CHROM")
            t = line.split('\t'); t.delete_at(t.size - 1); t[-1] = sample
            puts "##FILTER=<ID=HET1,Description=\"Heterozygous in the first haplotype\">"
            puts "##FILTER=<ID=HET2,Description=\"Heterozygous in the second haplotype\">"
            puts "##FILTER=<ID=GAP1,Description=\"Uncalled in the first haplotype\">"
            puts "##FILTER=<ID=GAP2,Description=\"Uncalled in the second haplotype\">"
            line = t.join('\t')
          end
          puts line; next
        end
        t = line.split('\t')
        next unless re_ctg.match(t[0])
        ad : Array(Int32)? = nil
        filters = [] of String; ht = [nil, nil] of String?
        [0, 1].each do |hap_i|
          m = /^(\.|[0-9]+)\/(\.|[0-9]+):(\S+)/.match(t[9 + hap_i])
          unless m
            STDERR.puts "Malformatted VCF"; return 1
          end
          s = m[3].split(',')
          the_ad = (ad ||= Array.new(s.size, 0))
          s.each_with_index { |v, j| the_ad[j] += v.to_i }
          if m[1] == '.'
            filters << "GAP#{hap_i + 1}"; ht[hap_i] = "."
          elsif m[1] != m[2]
            filters << "HET#{hap_i + 1}"; ht[hap_i] = "."
          else
            ht[hap_i] = m[1]
          end
        end
        t.delete_at(t.size - 1)
        hap = 0; st2 = t[1].to_i; en2 = st2 + t[3].size
        if is_male
          if t[0] =~ /^(chr)?X/
            if (hv = hgver) && (pars = par[hv]?)
              in_par = pars.any? { |par_r| par_r[0] <= st2 && en2 <= par_r[1] }
              hap = in_par ? 0 : 2
            end
          elsif t[0] =~ /^(chr)?Y/
            hap = 1
          end
        end
        if hap > 0 && filters.size == 1
          filters.clear if (hap == 2 && filters[0] == "GAP1") || (hap == 1 && filters[0] == "GAP2")
        end
        t[5] = "30"; t[6] = filters.empty? ? "." : filters.join(";")
        t << ht.map { |hap_val| hap_val || "." }.join("|") + ":" + (ad || [] of Int32).join(",")
        puts t.join('\t')
      end
    end
    0
  end
end
