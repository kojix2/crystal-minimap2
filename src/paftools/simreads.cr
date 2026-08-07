module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # mason2fq — convert mason2-simulated SAM to FASTQ
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  # Mirrors paftools.js's own `Bytes.prototype.revcomp` table (paftools.js
  # lines ~213-224) — NOT the unrelated global k8_revcomp C++ builtin. Only
  # characters in this exact IUPAC set are complemented; any other byte maps
  # to NUL (paftools.js initializes the 256-entry lookup table to 0 and only
  # overwrites entries for these 32 characters), so we replicate that default
  # instead of leaving unmapped characters unchanged.
  REVCOMP_TABLE = begin
    s1 = "WSATUGCYRKMBDHVNwsatugcyrkmbdhvn"
    s2 = "WSTAACGRYMKVHDBNwstaacgrymkvhdbn"
    h = Hash(Char, Char).new('\u0000')
    s1.each_char.with_index { |c, i| h[c] = s2[i] }
    h
  end

  private def self.revcomp_seq(s : String) : String
    String.build do |sb|
      (s.size - 1).downto(0) { |i| sb << REVCOMP_TABLE[s[i]] }
    end
  end

  # Fields for a pending mason2 SAM record, mirroring the JS `last` array:
  # [qname, chr, pos, pos_end, strand, seq, qual, read_no, comment]
  private struct Mason2Rec
    property qname : String
    property chr : String
    property pos : Int32
    property pos_end : Int32
    property strand : String
    property seq : String
    property qual : String
    property read_no : Int32
    property comment : String

    def initialize(@qname, @chr, @pos, @pos_end, @strand, @seq, @qual, @read_no, @comment)
    end
  end

  def self.cmd_mason2fq(args : Array(String)) : Int32
    if args.empty?
      STDERR.puts "Usage: paftools mason2fq <mason.sam>"
      return 1
    end

    print_se = ->(a : Mason2Rec) {
      puts "@#{a.qname}!#{a.chr}!#{a.pos}!#{a.pos_end}!#{a.strand} #{a.comment}"
      puts a.seq
      puts "+"
      puts a.qual
    }

    re_cigar = /(\d+)([MIDSHN])/
    last : Mason2Rec? = nil

    open_in(args[0]) do |io|
      io.each_line(chomp: true) do |line|
        next if line.starts_with?('@')
        t = line.split('\t')
        l_ref = 0
        t[5].scan(re_cigar) do |m|
          l_ref += m[1].to_i if m[2] == "D" || m[2] == "M" || m[2] == "N"
        end
        flag = t[1].to_i
        rev = (flag & 16) != 0
        seq : String; qual : String
        if rev
          seq = revcomp_seq(t[9])
          qual = t[10].reverse
        else
          seq = t[9]; qual = t[10]
        end
        # NOTE (bug fix): the JS original uses /^simulated./ — an unescaped
        # dot matches ANY character, so it would strip "simulatedX" prefixes
        # of any 10-char form, not just the literal string "simulated.". This
        # is almost certainly meant to be a literal dot; we use the escaped
        # regex here.
        qname = t[0].sub(/^simulated\./, "")
        chr = t[2]
        pos = t[3].to_i - 1
        strand = (flag & 16) != 0 ? "-" : "+"
        read_no_bits = flag & 0xc0
        read_no = read_no_bits == 0x40 ? 1 : read_no_bits == 0x80 ? 2 : 0
        err = 0; snp = 0; indel = 0
        (11...t.size).each do |j|
          if m = /^XE:i:(\d+)/.match(t[j])
            err = m[1].to_i
          elsif m = /^XS:i:(\d+)/.match(t[j])
            snp = m[1].to_i
          elsif m = /^XI:i:(\d+)/.match(t[j])
            indel = m[1].to_i
          end
        end
        comment = "#{err}:#{snp}:#{indel}"

        cur_last = last
        if cur_last.nil?
          last = Mason2Rec.new(qname, chr, pos, pos + l_ref, strand, seq, qual, read_no, comment)
        elsif cur_last.qname != qname
          print_se.call(cur_last)
          last = Mason2Rec.new(qname, chr, pos, pos + l_ref, strand, seq, qual, read_no, comment)
        else
          if read_no == 2 # cur_last is the first read
            raise "ERROR: can't find read1" if cur_last.read_no != 1
            name = "#{qname}!#{chr}!#{cur_last.pos}_#{pos}!#{cur_last.pos_end}_#{pos + l_ref}!#{cur_last.strand}#{strand}"
            puts "@#{name}/1 #{cur_last.comment}"; puts cur_last.seq; puts "+"; puts cur_last.qual
            puts "@#{name}/2 #{comment}"; puts seq; puts "+"; puts qual
          else
            raise "ERROR: can't find read2" if cur_last.read_no != 2
            name = "#{qname}!#{chr}!#{pos}_#{cur_last.pos}!#{pos + l_ref}_#{cur_last.pos_end}!#{strand}#{cur_last.strand}"
            puts "@#{name}/1 #{comment}"; puts seq; puts "+"; puts qual
            puts "@#{name}/2 #{cur_last.comment}"; puts cur_last.seq; puts "+"; puts cur_last.qual
          end
          last = nil
        end
      end
    end
    print_se.call(last.not_nil!) if last
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # sim2bed — convert simulated read names (mason2/pbsim/badread) to BED6
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_sim2bed(args : Array(String)) : Int32
    if args.empty?
      STDERR.puts "Usage: paftools sim2bed <sim.txt>"
      return 1
    end
    open_in(args[0]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('!')
        next if t.size < 5
        chr = t[1]
        st : Int32; en : Int32; strand : String
        if t[2].includes?('_') # mason paired-end
          pos = t[2].split('_')
          fin = t[3].split('_')
          m = /^(.)(.)\/([12])$/.match(t[4])
          next unless m
          strand = m[3] == "1" ? m[1] : m[2]
          read_no = m[3].to_i - 1
          st = pos[read_no].to_i
          en = fin[read_no].to_i
        else # badread/pbsim long reads
          st = t[2].to_i
          en = t[3].to_i
          strand = t[4]
        end
        if st > en
          st, en = en, st
        end
        puts [chr, st, en, line, 0, strand].join('\t')
      end
    end
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # pbsim2fq — convert pbsim MAF alignments to FASTA (name mirrors upstream;
  # the original also emits FASTA, not FASTQ, despite the command name)
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_pbsim2fq(args : Array(String)) : Int32
    if args.size < 2
      STDERR.puts "Usage: paftools pbsim2fq <ref.fa.fai> <pbsim1.maf> [[pbsim2.maf] ...]"
      return 1
    end

    chr_list = [] of String
    open_in(args[0]) do |io|
      io.each_line(chomp: true) do |line|
        chr_list << line.split(/\s+/)[0]
      end
    end

    (1...args.size).each do |k|
      fn = args[k]
      state = 0
      reg_st = 0; reg_en = 0
      open_in(fn) do |io|
        io.each_line(chomp: true) do |line|
          if state == 0 && line[0]? == 'a'
            state = 1
          elsif state == 1 && line[0]? == 's'
            t = line.split(/\s+/)
            st = t[2].to_i
            reg_st = st; reg_en = st + t[3].to_i
            state = 2
          elsif state == 2 && line[0]? == 's'
            t = line.split(/\s+/)
            m = /S(\d+)_\d+/.match(t[1])
            raise "Failed to parse the read name" unless m
            chr_id = m[1].to_i - 1
            raise "Index outside the chr list" if chr_id >= chr_list.size
            name = "#{t[1]}!#{chr_list[chr_id]}!#{reg_st}!#{reg_en}!#{t[4]}"
            seq = t[6].delete('-')
            raise "Inconsistent read length" if seq.size != t[5].to_i
            unless seq.includes?("NN")
              seq = revcomp_seq(seq) if t[4] == "-"
              puts ">#{name}"
              puts seq
            end
            state = 0
          end
        end
      end
    end
    0
  end

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # badread2fa — convert Badread-simulated FASTQ/FASTA headers to FASTA
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  def self.cmd_badread2fa(args : Array(String)) : Int32
    if args.size < 2
      STDERR.puts "Usage: paftools badread2fa <ref.fa.fai> <badread.fq>"
      return 1
    end

    len = Hash(String, Int32).new
    open_in(args[0]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('\t')
        len[t[0]] = t[1].to_i
      end
    end

    id = 0; n_discard = 0
    open_in(args[1]) do |io|
      loop do
        line = io.gets(chomp: true)
        break unless line
        is_fq = line[0]? == '@'
        tag = ""; a : Array(String)? = nil
        if !/\schimera\s/.matches?(line) &&
           (m = /\s(\S+),([+-])strand,(\d+)-(\d+).*read_identity=([0-9.]+)%/.match(line))
          contig = m[1]
          raise "failed to find the contig length of #{contig}" unless len.has_key?(contig)
          st = m[3].to_i; en = m[4].to_i
          a = if m[2] == "+"
                ["S#{id + 1}", contig, st.to_s, en.to_s, m[2]]
              else
                clen = len[contig]
                ["S#{id + 1}", contig, (clen - en).to_s, (clen - st).to_s, m[2]]
              end
          tag = "ri:f:#{m[5]}"
        end
        seq = io.gets(chomp: true) || ""
        if is_fq
          io.gets(chomp: true) # '+' line
          io.gets(chomp: true) # qual line
        end
        if av = a
          puts ">#{av.join("!")}\t#{tag}"
          puts seq
        else
          n_discard += 1
        end
        id += 1
      end
    end
    STDERR.puts "WARNING: discarded #{n_discard} reads"
    0
  end
end
