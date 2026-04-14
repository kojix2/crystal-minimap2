module Minimap2
  # ---------------------------------------------------------------------------
  # Nucleotide encoding table: A=0, C=1, G=2, T=3, other=4
  # ---------------------------------------------------------------------------
  SEQ_NT4_TABLE = begin
    t = Array(UInt8).new(256, 4_u8)
    # uppercase
    t['A'.ord] = 0_u8; t['C'.ord] = 1_u8
    t['G'.ord] = 2_u8; t['T'.ord] = 3_u8; t['U'.ord] = 3_u8
    # lowercase
    t['a'.ord] = 0_u8; t['c'.ord] = 1_u8
    t['g'.ord] = 2_u8; t['t'.ord] = 3_u8; t['u'.ord] = 3_u8
    t
  end

  # Complement table for seq_comp_table (used when reversing sequences)
  SEQ_COMP_TABLE = begin
    t = Array(UInt8).new(256, &.to_u8)
    pairs = {
      'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C',
      'a' => 't', 't' => 'a', 'c' => 'g', 'g' => 'c',
      'U' => 'A', 'u' => 'a',
      'R' => 'Y', 'Y' => 'R', 'S' => 'S', 'W' => 'W',
      'K' => 'M', 'M' => 'K', 'B' => 'V', 'V' => 'B',
      'D' => 'H', 'H' => 'D', 'N' => 'N',
      'r' => 'y', 'y' => 'r', 's' => 's', 'w' => 'w',
      'k' => 'm', 'm' => 'k', 'b' => 'v', 'v' => 'b',
      'd' => 'h', 'h' => 'd', 'n' => 'n',
    }
    pairs.each { |k, v| t[k.ord] = v.ord.to_u8 }
    t
  end

  # ---------------------------------------------------------------------------
  # Tiny circular queue  (mirrors tiny_queue_t from sketch.c)
  # Capacity is fixed at 32 (1 << 5).
  # ---------------------------------------------------------------------------
  private struct TinyQueue
    @front : Int32
    @count : Int32
    @a : StaticArray(Int32, 32)

    def initialize
      @front = 0
      @count = 0
      @a = StaticArray(Int32, 32).new(0)
    end

    def push(x : Int32)
      @a[(@count + @front) & 0x1f] = x
      @count += 1
    end

    def shift : Int32
      return -1 if @count == 0
      x = @a[@front]
      @front = (@front + 1) & 0x1f
      @count -= 1
      x
    end

    def size : Int32
      @count
    end

    def clear
      @front = @count = 0
    end
  end

  # ---------------------------------------------------------------------------
  # Invertible hash (from sketch.c — hash64)
  # ---------------------------------------------------------------------------
  @[AlwaysInline]
  private def self.hash64(key : UInt64, mask : UInt64) : UInt64
    key = (~key &+ (key << 21)) & mask
    key = key ^ (key >> 24)
    key = ((key &+ (key << 3)) &+ (key << 8)) & mask
    key = key ^ (key >> 14)
    key = ((key &+ (key << 2)) &+ (key << 4)) & mask
    key = key ^ (key >> 28)
    key = (key &+ (key << 31)) & mask
    key
  end

  # ---------------------------------------------------------------------------
  # Compute (w,k)-minimizers on a DNA sequence, appending to *p*.
  #
  # Encoding:
  #   p[i].x = hash<<8 | kmer_span
  #   p[i].y = rid<<32 | last_pos<<1 | strand
  #
  # Equivalent to mm_sketch() in sketch.c.
  # ---------------------------------------------------------------------------
  def self.mm_sketch(seq : String | Bytes, len : Int32, w : Int32, k : Int32,
                     rid : UInt32, is_hpc : Bool, p : Array(Mm128)) : Nil
    raise ArgumentError.new("len must be >0") unless len > 0
    raise ArgumentError.new("w must be in [1,255]") unless w > 0 && w < 256
    raise ArgumentError.new("k must be in [1,28]") unless k > 0 && k <= 28

    shift1 = 2 * (k - 1)
    mask = (1_u64 << (2 * k)) &- 1_u64
    kmer = StaticArray(UInt64, 2).new(0_u64)
    buf = Array(Mm128).new(w + 1) { Mm128.max } # circular buffer of size w
    min = Mm128.max
    tq = TinyQueue.new

    buf_pos = 0
    min_pos = 0
    kmer_span = 0
    l = 0

    seq_bytes = seq.is_a?(String) ? seq.to_slice : seq

    i = 0
    while i < len
      c = SEQ_NT4_TABLE[seq_bytes[i].to_i]
      info = Mm128.max

      if c < 4 # unambiguous base
        if is_hpc
          # skip homopolymer run
          skip_len = 1
          if i + 1 < len && SEQ_NT4_TABLE[seq_bytes[i + 1].to_i] == c
            skip_len = 2
            while i + skip_len < len && SEQ_NT4_TABLE[seq_bytes[i + skip_len].to_i] == c
              skip_len += 1
            end
          end
          i += skip_len - 1 # advance to end of run
          tq.push(skip_len)
          kmer_span += skip_len
          kmer_span -= tq.shift if tq.size > k
        else
          kmer_span = l + 1 < k ? l + 1 : k
        end

        kmer[0] = ((kmer[0] << 2) | c.to_u64) & mask              # forward k-mer
        kmer[1] = (kmer[1] >> 2) | ((3_u64 ^ c.to_u64) << shift1) # reverse k-mer
        # skip palindromic k-mers — strand unknown
        unless kmer[0] == kmer[1]
          z = kmer[0] < kmer[1] ? 0 : 1 # strand
          l += 1
          if l >= k && kmer_span < 256
            info = Mm128.new(
              hash64(kmer[z], mask) << 8 | kmer_span.to_u64,
              rid.to_u64 << 32 | (i.to_u32).to_u64 << 1 | z.to_u64
            )
          end
        end
      else
        # ambiguous base: reset
        l = 0
        tq.clear
        kmer_span = 0
      end

      buf[buf_pos] = info

      # First window: flush identical k-mers that were deferred
      if l == w + k - 1 && !min.x.==(Mm128::UINT64_MAX)
        (buf_pos + 1).upto(w - 1) do |j|
          p << buf[j] if buf[j].x == min.x && buf[j].y != min.y
        end
        0.upto(buf_pos - 1) do |j|
          p << buf[j] if buf[j].x == min.x && buf[j].y != min.y
        end
      end

      if info.x <= min.x
        # new minimum — emit old minimum first
        p << min if l >= w + k && min.x != Mm128::UINT64_MAX
        min = info
        min_pos = buf_pos
      elsif buf_pos == min_pos
        # old minimum slid out of window — find new minimum
        p << min if l >= w + k - 1 && min.x != Mm128::UINT64_MAX
        # scan entire window for new minimum
        min = Mm128.max
        min_pos = 0
        (buf_pos + 1).upto(w - 1) do |j|
          if min.x >= buf[j].x
            min = buf[j]
            min_pos = j
          end
        end
        0.upto(buf_pos) do |j|
          if min.x >= buf[j].x
            min = buf[j]
            min_pos = j
          end
        end
        # emit identical k-mers to the new minimum
        if l >= w + k - 1 && min.x != Mm128::UINT64_MAX
          (buf_pos + 1).upto(w - 1) do |j|
            p << buf[j] if buf[j].x == min.x && min.y != buf[j].y
          end
          0.upto(buf_pos) do |j|
            p << buf[j] if buf[j].x == min.x && min.y != buf[j].y
          end
        end
      end

      buf_pos += 1
      buf_pos = 0 if buf_pos == w

      i += 1
    end

    p << min if min.x != Mm128::UINT64_MAX
  end

  # Overload accepting Slice(UInt8) directly (internal use)
  def self.mm_sketch(seq : Slice(UInt8), w : Int32, k : Int32,
                     rid : UInt32, is_hpc : Bool, p : Array(Mm128)) : Nil
    mm_sketch(seq, seq.size, w, k, rid, is_hpc, p)
  end
end
