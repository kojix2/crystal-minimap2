module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # gff2bed — convert GTF/GFF3 to BED12 (or junction BED with -j)
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_gff2bed(args : Array(String)) : Int32
    is_short = false; print_junc = false; output_gene = false; ens_canon_only = false
    fn_ucsc : String? = nil
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-s"          ; is_short = true
      when "-j"          ; print_junc = true
      when "-G"          ; output_gene = true
      when "-e"          ; ens_canon_only = true
      when "-u"          ; i += 1; fn_ucsc = args[i]
      when "-h", "--help"; STDERR.puts "Usage: paftools gff2bed [-sjGe] [-u fai] <in.gff>"; return 0
      else                 rest << args[i]
      end
      i += 1
    end
    if rest.empty?
      STDERR.puts "Usage: paftools gff2bed [-sjGe] [-u fai] <in.gff>"; return 1
    end

    ens2ucsc = Hash(String, String).new
    if fai = fn_ucsc
      File.open(fai) do |file|
        file.each_line(chomp: true) do |line|
          t = line.split('\t'); s = t[0]
          if s =~ /_(random|alt|decoy)$/
            s = s.gsub(/_(random|alt|decoy)$/, "").gsub(/^chr\S+_/, "")
          else
            s = s.gsub(/^chrUn_/, "")
          end
          s = s.gsub(/v(\d+)/, ".\\1")
          ens2ucsc[s] = t[0] if s != t[0]
        end
      end
    end

    colors = {"protein_coding" => "0,128,255", "mRNA" => "0,128,255",
              "lincRNA" => "0,192,0", "snRNA" => "0,192,0",
              "miRNA" => "0,192,0", "misc_RNA" => "0,192,0"}

    print_bed12 = ->(exons : Array(Array(String)), cds_st : Int32, cds_en : Int32) {
      return if exons.empty?
      name = is_short ? "#{exons[0][7]}|#{exons[0][5]}" : exons[0][4..6].join("|")
      sorted = exons.sort_by(&.[1].to_i)
      if print_junc
        (1...sorted.size).each do |j|
          puts [sorted[j][0], sorted[j - 1][2], sorted[j][1], name, 1000, sorted[j][3]].join('\t')
        end
        return
      end
      st = sorted[0][1].to_i; en = sorted.last[2].to_i
      cds_st = st if cds_st == (1 << 30); cds_en = en if cds_en == 0
      sizes = sorted.map { |e| e[2].to_i - e[1].to_i }
      starts = sorted.map { |e| e[1].to_i - st }
      color = colors[sorted[0][5]]? || "196,196,196"
      puts [sorted[0][0], st, en, name, 1000, sorted[0][3], cds_st, cds_en, color,
            sorted.size, sizes.join(",") + ",", starts.join(",") + ","].join('\t')
    }

    re_gtf = /\b(transcript_id|transcript_type|transcript_biotype|gene_name|gene_id|gbkey|transcript_name|tag) "([^"]+)";/
    re_gff3 = /\b(transcript_id|transcript_type|transcript_biotype|gene_name|gene_id|gbkey|transcript_name)=([^;]+)/
    re_gtf_gene = /\b(gene_id|gene_type|gene_name) "([^;]+)";/
    re_gff3_gene = /\b(gene_id|gene_type|source_gene|gene_biotype|gene_name)=([^;]+);/

    exons = [] of Array(String); cds_st = 1 << 30; cds_en = 0; last_id : String? = nil

    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('\t')
        next if t.size < 9 || t[0].starts_with?('#')

        if fn_ucsc
          t[0] = ens2ucsc[t[0]]? || t[0]
          puts t.join('\t'); next
        end

        if output_gene
          next if t[2] != "gene"
          id = nil; type = ""; name = "N/A"
          line.scan(re_gtf_gene) { |mat| id = mat[2] if mat[1] == "gene_id"; type = mat[2] if mat[1] == "gene_type"; name = mat[2] if mat[1] == "gene_name" }
          line.scan(re_gff3_gene) { |mat| id = mat[2] if mat[1] == "gene_id" || mat[1] == "source_gene"; type = mat[2] if mat[1] == "gene_type" || mat[1] == "gene_biotype"; name = mat[2] if mat[1] == "gene_name" }
          puts [t[0], t[3].to_i - 1, t[4], "#{id}|#{type}|#{name}", 1000, t[6]].join('\t')
          next
        end

        next if t[2] != "CDS" && t[2] != "exon"
        t[3] = (t[3].to_i - 1).to_s; t[4] = t[4].to_i.to_s

        id = nil; type = ""; name = "N/A"; biotype = ""
        tname = "N/A"; ens_canonical = false
        line.scan(re_gtf) do |mat|
          case mat[1]
          when "transcript_id"              ; id = mat[2]
          when "transcript_type"            ; type = mat[2]
          when "transcript_biotype", "gbkey"; biotype = mat[2]
          when "gene_name", "gene_id"       ; name = mat[2]
          when "transcript_name"            ; tname = mat[2]
          when "tag"                        ; ens_canonical = true if mat[2] == "Ensembl_canonical"
          end
        end
        line.scan(re_gff3) do |mat|
          case mat[1]
          when "transcript_id"              ; id = mat[2]
          when "transcript_type"            ; type = mat[2]
          when "transcript_biotype", "gbkey"; biotype = mat[2]
          when "gene_name", "gene_id"       ; name = mat[2]
          when "transcript_name"            ; tname = mat[2]
          end
        end
        next if ens_canon_only && !ens_canonical
        type = biotype if type.empty? && !biotype.empty?
        next unless id
        if id != last_id
          print_bed12.call(exons, cds_st, cds_en)
          exons = [] of Array(String); cds_st = 1 << 30; cds_en = 0; last_id = id
        end
        chr = ens2ucsc[t[0]]? || t[0]
        chr = chr.gsub(/([A-Z]+\d+)\.(\d+)/, "chrUn_\\1v\\2") if fn_ucsc && chr =~ /^[A-Z]+\d+\.\d+$/
        if t[2] == "CDS"
          cds_st = [cds_st, t[3].to_i].min; cds_en = [cds_en, t[4].to_i].max
        else
          exons << [chr, t[3], t[4], t[6], id || "", type, name, tname]
        end
      end
    end
    print_bed12.call(exons, cds_st, cds_en) if last_id
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # gff2junc — extract splice junctions from GFF3
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_gff2junc(args : Array(String)) : Int32
    feat = "CDS"; rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-f"          ; i += 1; feat = args[i]
      when "-h", "--help"; STDERR.puts "Usage: paftools gff2junc [-f feature] <in.gff3>"; return 0
      else                 rest << args[i]
      end
      i += 1
    end
    if rest.empty?
      STDERR.puts "Usage: paftools gff2junc [-f feature] <in.gff3>"; return 1
    end

    process = ->(a : Array(Array(String))) {
      return if a.size < 2
      s = a.sort_by(&.[4].to_i)
      (1...s.size).each { |j| puts [s[j][1], s[j - 1][5], s[j][4], s[j][0], 0, s[j][7]].join('\t') }
    }

    buf = [] of Array(String)
    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('\t'); next if t.size < 9 || t[0][0] == '#'
        next if t[2].downcase != feat.downcase
        m = /\bParent=([^;]+)/.match(t[8])
        unless m
          STDERR.puts "Can't find Parent"; next
        end
        t[3] = (t[3].to_i - 1).to_s
        entry = [m[1], t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8]]
        if !buf.empty? && buf[0][0] != m[1]
          process.call(buf); buf.clear
        end
        buf << entry
      end
    end
    process.call(buf)
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # splice2bed — splice alignments (PAF/SAM) → BED12
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_splice2bed(args : Array(String)) : Int32
    keep_multi = false; fn_conv : String? = nil
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-m"          ; keep_multi = true
      when "-n"          ; i += 1; fn_conv = args[i]
      when "-h", "--help"; STDERR.puts "Usage: paftools splice2bed [-m] <in.paf>|<in.sam>"; return 0
      else                 rest << args[i]
      end
      i += 1
    end
    if rest.empty?
      STDERR.puts "Usage: paftools splice2bed [-m] <in.paf>|<in.sam>"; return 1
    end

    conv = if fn = fn_conv
             Hash(String, String).new.tap do |hsh|
               File.open(fn) { |file| file.each_line(chomp: true) { |line2| t = line2.split('\t'); hsh[t[0]] = t[1] } }
             end
           end

    colors = ["0,128,255", "255,0,0", "0,192,0"]
    print_lines = ->(a : Array(Array(String))) {
      return if a.empty?
      n_pri = a.count { |row| row[8] == "0" }
      a.each { |row| row[8] = "1" } if n_pri > 1
      STDERR.puts "Warning: #{a[0][3]} has no primary alignment" if n_pri == 0
      a.each do |row|
        next if !keep_multi && row[8] == "2"
        row[8] = colors[row[8].to_i]
        puts row.join('\t')
      end
    }

    buf = [] of Array(String)
    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        next if line.starts_with?('@')
        t = line.split('\t')
        qname = conv.try { |cov_map| cov_map[t[0]]? } || t[0]
        t[0] = qname

        is_pri = false; cigar : String? = nil; a1 : Array(String)? = nil

        if t.size >= 12 && (t[4] == "+" || t[4] == "-") # PAF
          t[12..].each do |tag|
            if tag.starts_with?("cg:Z:")
              cigar = tag[5..]
            elsif tag.starts_with?("s2:i:")
              is_pri = true
            end
          end
          a1 = [t[5], t[7], t[8], t[0], (t[9].to_f / t[10].to_f * 1000).to_i.to_s, t[4]]
        elsif t.size >= 10 # SAM
          flag = t[1].to_i
          next if (flag & 4) != 0
          cigar = t[5]; is_pri = (flag & 0x100) == 0
          if (flag & 1) != 0
            t[0] += "/#{(flag >> 6) & 3}"
          end
          te_val : String? = nil
          a1 = [t[2], (t[3].to_i - 1).to_s, te_val || "?", t[0], "1000", (flag & 16) != 0 ? "-" : "+"]
        else
          next
        end

        if !buf.empty? && buf[0][3] != t[0]
          print_lines.call(buf); buf.clear
        end

        next unless cigar
        x0 = 0; x = 0; bs = [] of Int32; bl = [] of Int32
        cigar.scan(/(\d+)([MIDNSHP=X])/) do |mat|
          l = mat[1].to_i; op = mat[2][0]
          case op
          when 'M', 'D'; x += l
          when 'N'
            bs << x0; bl << (x - x0); x0 = x + l; x += l
          when 'S', 'H' # ignore for BED
          end
        end
        bs << x0; bl << (x - x0)

        row = (a1 || next).dup
        if row[2] == "?"
          row[2] = (row[1].to_i + x).to_s
        end
        row.concat([row[1], row[2], is_pri ? "0" : "2",
                    bs.size.to_s, bl.map(&.to_s).join(",") + ",",
                    bs.map(&.to_s).join(",") + ","])
        buf << row
      end
    end
    print_lines.call(buf)
    0
  end
end
