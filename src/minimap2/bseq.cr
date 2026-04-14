require "compress/gzip"

module Minimap2
  # ---------------------------------------------------------------------------
  # Biological sequence record — mirrors mm_bseq1_t in bseq.h
  # ---------------------------------------------------------------------------
  class BSeq1
    property name : String
    property seq : String
    property qual : String?
    property comment : String?
    property l_seq : Int32
    property rid : UInt32 # reference id assigned during index construction

    def initialize(@name, @seq, @qual = nil, @comment = nil)
      @l_seq = @seq.size
      @rid = 0_u32
    end

    # Convert U/u → T/t in-place (FASTA convention)
    def normalize_u_to_t! : Nil
      @seq = @seq.gsub('U', 'T').gsub('u', 't')
    end
  end

  # ---------------------------------------------------------------------------
  # FASTA/FASTQ file reader — pure Crystal replacement for bseq.c / kseq.h.
  # Supports plain text and gzip-compressed files (via Crystal's Compress::Gzip).
  # ---------------------------------------------------------------------------
  class BSeqFile
    CHECK_PAIR_THRES = 1_000_000

    @io : IO
    @pending : BSeq1? # stashed record that overflowed the last chunk

    def self.open(fn : String) : BSeqFile
      io : IO
      if fn == "-"
        io = STDIN
      elsif fn.ends_with?(".gz") || fn.ends_with?(".bgz")
        io = Compress::Gzip::Reader.new(File.open(fn))
      else
        io = File.open(fn)
      end
      new(io)
    end

    def initialize(@io : IO)
      @pending = nil
    end

    def close
      @io.close
    end

    def eof? : Bool
      @pending.nil? && @io.peek.try(&.empty?) != false
    end

    # Read a mini-batch of sequences up to *chunk_size* total bases.
    # Returns an empty array at EOF.
    # *with_qual*    — keep quality string (FASTQ only)
    # *with_comment* — keep FASTQ/FASTA comment field
    # *frag_mode*    — keep paired reads together even if chunk overflows
    def read_seqs(chunk_size : Int64, with_qual : Bool = false,
                  with_comment : Bool = false, frag_mode : Bool = false) : Array(BSeq1)
      a = [] of BSeq1
      size = 0_i64

      if p = @pending
        a << p
        size += p.l_seq
        @pending = nil
      end

      while s = next_seq(with_qual, with_comment)
        a << s
        size += s.l_seq
        if size >= chunk_size
          if frag_mode && a.last.l_seq < CHECK_PAIR_THRES
            # keep reading until we get a read with a different name (pair boundary)
            while s2 = next_seq(with_qual, with_comment)
              if same_qname?(s2.name, a.last.name)
                a << s2
              else
                @pending = s2
                break
              end
            end
          end
          break
        end
      end

      a
    end

    # ---------------------------------------------------------------------------
    # Parse one FASTA or FASTQ record from @io.
    # Returns nil at end of file.
    # ---------------------------------------------------------------------------
    private def next_seq(with_qual : Bool, with_comment : Bool) : BSeq1?
      # Skip blank lines and lines that don't start a record
      line = read_non_empty_line
      return if line.nil?

      if line.starts_with?('>')
        # FASTA
        header = line[1..]
        name, comment = split_header(header, with_comment)
        seq = read_fasta_sequence
        s = BSeq1.new(name, seq, nil, with_comment ? comment : nil)
        s.normalize_u_to_t!
        STDERR.puts "[WARNING] empty sequence name in the input." if name.empty?
        s
      elsif line.starts_with?('@')
        # FASTQ
        header = line[1..]
        name, comment = split_header(header, with_comment)
        seq_line = @io.gets(chomp: true) || ""
        seq = seq_line.strip
        @io.gets(chomp: true) # consume '+' line
        qual_line = @io.gets(chomp: true) || ""
        qual = qual_line.strip
        STDERR.puts "[WARNING] empty sequence name in the input." if name.empty?
        s = BSeq1.new(
          name, seq.gsub('U', 'T').gsub('u', 't'),
          (with_qual && !qual.empty?) ? qual : nil,
          with_comment ? comment : nil
        )
        s
      else
        # Unexpected format — skip line and try again
        next_seq(with_qual, with_comment)
      end
    end

    # Read a FASTA sequence body (possibly multi-line) until next '>' or EOF
    private def read_fasta_sequence : String
      buf = String::Builder.new
      loop do
        peeked = @io.peek
        break if peeked.nil? || peeked.empty?
        break if peeked[0] == '>'.ord || peeked[0] == '@'.ord
        line = @io.gets(chomp: true)
        break if line.nil?
        stripped = line.strip
        buf << stripped unless stripped.empty?
      end
      buf.to_s
    end

    # Skip blank lines; return first non-empty line or nil
    private def read_non_empty_line : String?
      loop do
        line = @io.gets(chomp: true)
        return if line.nil?
        stripped = line.strip
        return stripped unless stripped.empty?
      end
    end

    # Split "name [comment]" on first whitespace
    private def split_header(header : String, want_comment : Bool) : {String, String?}
      i = header.each_char.with_index.find { |chr, _| chr == ' ' || chr == '\t' }.try(&.[1])
      if i
        name = header[0...i]
        comment = want_comment ? header[(i + 1)..].lstrip : nil
        {name, comment}
      else
        {header.strip, nil}
      end
    end

    # Check if two query names differ only in /1 /2 suffix (bseq.c mm_qname_same)
    private def same_qname?(a : String, b : String) : Bool
      return false if a.empty? || b.empty?
      la = a.size; lb = b.size
      # strip trailing /1 or /2
      ea = (la >= 2 && a[la - 2] == '/' && (a[la - 1] == '1' || a[la - 1] == '2')) ? la - 2 : la
      eb = (lb >= 2 && b[lb - 2] == '/' && (b[lb - 1] == '1' || b[lb - 1] == '2')) ? lb - 2 : lb
      return false if ea != eb
      a[0...ea] == b[0...eb]
    end
  end

  # Multi-file fragment reader (mirrors mm_bseq_read_frag2)
  def self.bseq_read_frag(fps : Array(BSeqFile), chunk_size : Int64,
                          with_qual : Bool = false, with_comment : Bool = false) : Array(BSeq1)
    a = [] of BSeq1
    size = 0_i64

    loop do
      seqs = fps.map(&.send(:next_seq, with_qual, with_comment))
      if seqs.any?(Nil)
        if seqs.any? { |seq| !seq.nil? }
          STDERR.puts "[WARNING] query files have different number of records; extra records skipped."
        end
        break
      end
      seqs.each { |seq| next unless seq; a << seq; size += seq.l_seq }
      break if size >= chunk_size
    end
    a
  end
end
