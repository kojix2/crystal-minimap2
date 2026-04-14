module Minimap2
  # ---------------------------------------------------------------------------
  # One bucket of the hash index
  # During construction: a[] collects (hash, position) minimizer pairs.
  # After worker_post: a[] is cleared and h[] is the lookup table.
  # ---------------------------------------------------------------------------
  private class MmIdxBucket
    property a : Array(Mm128)                 # collecting phase
    property h : Hash(UInt64, Array(UInt64))  # minimizer_key → sorted positions

    def initialize
      @a = [] of Mm128
      @h = {} of UInt64 => Array(UInt64)
    end
  end

  # ---------------------------------------------------------------------------
  # The minimap2 index — mirrors mm_idx_t
  # ---------------------------------------------------------------------------
  class MmIdx
    property w      : Int32   # window size
    property k      : Int32   # k-mer size
    property b      : Int32   # bucket bits
    property flag   : Int32   # MM_I_* flags
    property n_seq  : UInt32
    property seq    : Array(IdxSeq)
    property s_arr  : Array(UInt32)   # packed 4-bit sequence (S in C)
    property name_map : Hash(String, UInt32)  # name → index
    property index  : Int32           # part index (for split index)
    property n_alt  : Int32           # number of ALT sequences
    @b_arr : Array(MmIdxBucket)

    protected getter b_arr : Array(MmIdxBucket)

    def initialize(@w, @k, @b, @flag)
      @b = [@b, @k * 2].min
      @w = [@w, 1].max
      @n_seq = 0_u32
      @seq = [] of IdxSeq
      @s_arr = [] of UInt32
      @name_map = {} of String => UInt32
      @index = 0
      @n_alt = 0
      @b_arr = Array(MmIdxBucket).new(1 << @b) { MmIdxBucket.new }
    end

    # -------------------------------------------------------------------------
    # Sequence packing helpers (mm_seq4_set / mm_seq4_get)
    # 8 bases per UInt32, 4 bits each.
    # -------------------------------------------------------------------------
    def seq4_set(offset : UInt64, c : Int32) : Nil
      idx  = (offset >> 3).to_i32
      @s_arr[idx] = @s_arr[idx] | (c.to_u32 << ((offset & 7) << 2))
    end

    def seq4_get(offset : UInt64) : UInt8
      ((@s_arr[(offset >> 3).to_i32] >> ((offset & 7) << 2)) & 0xf_u32).to_u8
    end

    # -------------------------------------------------------------------------
    # Ensure s_arr has room for total_len bases (allocating in chunks of 8).
    # -------------------------------------------------------------------------
    def ensure_seq_capacity(total_len : UInt64) : Nil
      needed = ((total_len + 7) / 8).to_i32
      if @s_arr.size < needed
        @s_arr.concat(Array(UInt32).new(needed - @s_arr.size, 0_u32))
      end
    end

    # -------------------------------------------------------------------------
    # Add minimizers to buckets (mirrors mm_idx_add).
    # -------------------------------------------------------------------------
    def add(a : Array(Mm128)) : Nil
      mask = (1 << @b) - 1
      a.each do |m|
        @b_arr[((m.x >> 8) & mask).to_i32].a << m
      end
    end

    # -------------------------------------------------------------------------
    # Sort and hash one bucket (mirrors worker_post).
    # -------------------------------------------------------------------------
    private def build_bucket(bi : Int32) : Nil
      bkt = @b_arr[bi]
      return if bkt.a.empty?

      # Sort by minimizer hash
      bkt.a.sort_by!(&.x)

      # Group by minimizer and build hash table
      h = bkt.h
      j = 0
      while j < bkt.a.size
        x0 = bkt.a[j].x >> 8
        # Find run of identical minimizers
        k = j + 1
        while k < bkt.a.size && bkt.a[k].x >> 8 == x0
          k += 1
        end
        key = x0 >> @b  # strip bucket-index bits
        positions = bkt.a[j...k].map(&.y).sort!
        h[key] = positions
        j = k
      end

      # Clear the construction array
      bkt.a.clear
    end

    # -------------------------------------------------------------------------
    # Process all buckets after collection (mirrors mm_idx_post).
    # Each bucket is independent, so we fan them out across a Parallel context.
    # -------------------------------------------------------------------------
    def post(n_threads : Int32 = 1) : Nil
      total = 1 << @b
      if n_threads <= 1
        total.times { |i| build_bucket(i) }
        return
      end

      ctx     = Fiber::ExecutionContext::Parallel.new("idx-post", n_threads)
      pending = Atomic(Int32).new(total)
      done_ch = Channel(Nil).new(1)

      total.times do |i|
        ctx.spawn do
          build_bucket(i)
          # sub returns the value *before* the decrement; when it was 1 we are last
          done_ch.send(nil) if pending.sub(1, :sequentially_consistent) == 1
        end
      end

      done_ch.receive
    end

    # -------------------------------------------------------------------------
    # Lookup positions for a minimizer (mirrors mm_idx_get).
    # Returns the array of hit positions, or nil if not found.
    # -------------------------------------------------------------------------
    def get(minier : UInt64) : Array(UInt64)?
      mask = ((1 << @b) - 1).to_u64
      key  = minier >> @b
      @b_arr[(minier & mask).to_i32].h[key]?
    end

    # -------------------------------------------------------------------------
    # Maximum occurrence threshold at fraction f (mirrors mm_idx_cal_max_occ).
    # -------------------------------------------------------------------------
    def cal_max_occ(f : Float32) : Int32
      return Int32::MAX if f <= 0.0_f32
      counts = [] of UInt32
      @b_arr.each do |bkt|
        bkt.h.each_value { |pos| counts << pos.size.to_u32 }
      end
      return Int32::MAX if counts.empty?
      (Minimap2.ks_ksmall_uint32(counts, ((1.0 - f) * counts.size).to_i32) + 1).to_i32
    end

    # -------------------------------------------------------------------------
    # Retrieve decoded sequence (mirrors mm_idx_getseq).
    # Returns encoded bases (0-3), or -1 on range error.
    # -------------------------------------------------------------------------
    def getseq(rid : UInt32, st : UInt32, en : UInt32, buf : Array(UInt8)) : Int32
      return -1 if rid >= @n_seq || st >= @seq[rid].len
      en = @seq[rid].len if en > @seq[rid].len
      off = @seq[rid].offset
      (st...en).each do |j|
        buf[j - st] = seq4_get(off + j)
      end
      (en - st).to_i32
    end

    def getseq_rev(rid : UInt32, st : UInt32, en : UInt32, buf : Array(UInt8)) : Int32
      return -1 if rid >= @n_seq || st >= @seq[rid].len
      s = @seq[rid]
      en = s.len if en > s.len
      st1 = s.offset + (s.len - en)
      en1 = s.offset + (s.len - st)
      (st1...en1).each do |i|
        c = seq4_get(i)
        buf[(en1 - i - 1).to_i32] = c < 4 ? (3_u8 - c) : c
      end
      (en - st).to_i32
    end

    def getseq2(is_rev : Bool, rid : UInt32, st : UInt32, en : UInt32, buf : Array(UInt8)) : Int32
      is_rev ? getseq_rev(rid, st, en, buf) : getseq(rid, st, en, buf)
    end

    # -------------------------------------------------------------------------
    # Build name → id map (mirrors mm_idx_index_name).
    # Returns true if duplicate names detected.
    # -------------------------------------------------------------------------
    def index_name : Bool
      has_dup = false
      @seq.each_with_index do |s, i|
        if @name_map.has_key?(s.name)
          has_dup = true
        else
          @name_map[s.name] = i.to_u32
        end
      end
      if has_dup && Minimap2.verbose >= 2
        STDERR.puts "[WARNING] some database sequences have identical sequence names"
      end
      has_dup
    end

    def name2id(name : String) : Int32
      v = @name_map[name]?
      v ? v.to_i32 : -1
    end

    # -------------------------------------------------------------------------
    # Print statistics (mirrors mm_idx_stat).
    # -------------------------------------------------------------------------
    def stat : Nil
      hpc = (@flag & I_HPC) != 0
      total_len = @seq.sum(&.len.to_u64)
      n_distinct = @b_arr.sum { |bkt| bkt.h.size.to_i64 }
      n_singletons = @b_arr.sum { |bkt| bkt.h.count { |_, v| v.size == 1 }.to_i64 }
      total_occ = @b_arr.sum { |bkt| bkt.h.sum { |_, v| v.size.to_i64 } }
      t = Minimap2.realtime - Minimap2.realtime0
      STDERR.printf(
        "[M::mm_idx_stat] kmer size: %d; skip: %d; is_hpc: %d; #seq: %d\n",
        @k, @w, hpc ? 1 : 0, @n_seq
      )
      if n_distinct > 0
        STDERR.printf(
          "[M::mm_idx_stat::%.3f*%.2f] distinct minimizers: %ld (%.2f%% singletons); avg occ: %.3f; avg spacing: %.3f; total len: %ld\n",
          t, Minimap2.cputime / t, n_distinct, 100.0 * n_singletons / n_distinct,
          total_occ.to_f / n_distinct, total_len.to_f / total_occ, total_len
        )
      end
    end

    # -------------------------------------------------------------------------
    # Serialise index to an IO (mirrors mm_idx_dump).
    # -------------------------------------------------------------------------
    def dump(io : IO) : Nil
      io.write(IDX_MAGIC.to_slice)
      sum_len = @seq.sum(&.len.to_u64)
      io.write_bytes(@w.to_u32, IO::ByteFormat::LittleEndian)
      io.write_bytes(@k.to_u32, IO::ByteFormat::LittleEndian)
      io.write_bytes(@b.to_u32, IO::ByteFormat::LittleEndian)
      io.write_bytes(@n_seq,    IO::ByteFormat::LittleEndian)
      io.write_bytes(@flag.to_u32, IO::ByteFormat::LittleEndian)
      @seq.each do |s|
        name_b = s.name.to_slice
        io.write_byte(name_b.size.to_u8)
        io.write(name_b)
        io.write_bytes(s.len, IO::ByteFormat::LittleEndian)
      end
      (1 << @b).times do |i|
        bkt = @b_arr[i]
        # Rebuild flat p[] and hash for serialisation
        flat_p  = [] of UInt64
        entries = [] of {UInt64, UInt64}
        bkt.h.each do |key, positions|
          if positions.size == 1
            # singleton: key (with singleton bit) → y_value
            entries << {(key << @b) << 1 | 1_u64, positions[0]}
          else
            start_p = flat_p.size.to_u64
            flat_p.concat(positions)
            entries << {(key << @b) << 1, (start_p << 32) | positions.size.to_u64}
          end
        end
        io.write_bytes(flat_p.size.to_u32, IO::ByteFormat::LittleEndian)
        flat_p.each { |v| io.write_bytes(v, IO::ByteFormat::LittleEndian) }
        io.write_bytes(entries.size.to_u32, IO::ByteFormat::LittleEndian)
        entries.each do |(k, v)|
          io.write_bytes(k, IO::ByteFormat::LittleEndian)
          io.write_bytes(v, IO::ByteFormat::LittleEndian)
        end
      end
      if (@flag & I_NO_SEQ) == 0
        n_words = ((sum_len + 7) / 8).to_i32
        n_words.times { |i| io.write_bytes(@s_arr[i]? || 0_u32, IO::ByteFormat::LittleEndian) }
      end
    end

    # -------------------------------------------------------------------------
    # Load index from IO (mirrors mm_idx_load).
    # -------------------------------------------------------------------------
    def self.load(io : IO) : MmIdx?
      magic = Bytes.new(4)
      return nil if io.read(magic) < 4
      return nil unless magic == IDX_MAGIC.to_slice
      w    = io.read_bytes(UInt32, IO::ByteFormat::LittleEndian).to_i32
      k    = io.read_bytes(UInt32, IO::ByteFormat::LittleEndian).to_i32
      b    = io.read_bytes(UInt32, IO::ByteFormat::LittleEndian).to_i32
      nseq = io.read_bytes(UInt32, IO::ByteFormat::LittleEndian)
      flag = io.read_bytes(UInt32, IO::ByteFormat::LittleEndian).to_i32
      mi = new(w, k, b, flag)
      mi.n_seq = nseq
      sum_len = 0_u64
      nseq.times do |_|
        l = io.read_byte || 0_u8
        name = l > 0 ? String.new(Bytes.new(l).tap { |buf| io.read_fully(buf) }) : ""
        len  = io.read_bytes(UInt32, IO::ByteFormat::LittleEndian)
        mi.seq << IdxSeq.new(name, sum_len, len)
        sum_len += len
      end
      (1 << b).times do |i|
        bkt = mi.b_arr[i]
        pn = io.read_bytes(UInt32, IO::ByteFormat::LittleEndian).to_i32
        flat_p = Array(UInt64).new(pn) { io.read_bytes(UInt64, IO::ByteFormat::LittleEndian) }
        size = io.read_bytes(UInt32, IO::ByteFormat::LittleEndian).to_i32
        size.times do
          raw_key = io.read_bytes(UInt64, IO::ByteFormat::LittleEndian)
          val     = io.read_bytes(UInt64, IO::ByteFormat::LittleEndian)
          key = raw_key >> b + 1  # strip bucket bits and singleton bit
          if (raw_key & 1) == 1
            bkt.h[key] = [val]
          else
            count   = (val & 0xffffffff_u64).to_i32
            start_p = (val >> 32).to_i64
            bkt.h[key] = flat_p[start_p...(start_p + count)]
          end
        end
      end
      if (flag & I_NO_SEQ) == 0
        n_words = ((sum_len + 7) / 8).to_i32
        mi.s_arr = Array(UInt32).new(n_words) { io.read_bytes(UInt32, IO::ByteFormat::LittleEndian) }
      end
      mi
    rescue IO::EOFError
      nil
    end

    # -------------------------------------------------------------------------
    # Check if a file is a prebuilt index; returns file size if yes, 0 if no.
    # -------------------------------------------------------------------------
    def self.is_idx?(fn : String) : Int64
      return 0_i64 if fn == "-"
      return -1_i64 unless File.exists?(fn)
      magic = Bytes.new(4)
      File.open(fn, "rb") do |f|
        return 0_i64 if f.read(magic) < 4
        return 0_i64 unless magic == IDX_MAGIC.to_slice[0...4]
        return f.size.to_i64
      end
    rescue
      -1_i64
    end

    # -------------------------------------------------------------------------
    # Build index from a BSeqFile pipeline (mirrors mm_idx_gen).
    # -------------------------------------------------------------------------
    def self.from_file(fp : BSeqFile, w : Int32, k : Int32, b : Int32,
                       flag : Int32, mini_batch_size : Int32,
                       n_threads : Int32, batch_size : UInt64) : MmIdx?
      return nil if fp.eof?

      mi = new(w, k, b, flag)
      sum_len = 0_u64
      p = [] of Mm128

      until fp.eof? || sum_len >= batch_size
        seqs = fp.read_seqs(mini_batch_size.to_i64)
        break if seqs.empty?

        # Assign RIDs and store metadata
        seqs.each do |s|
          seq_rec = IdxSeq.new(s.name, sum_len, s.l_seq.to_u32)
          mi.seq << seq_rec
          mi.n_seq += 1
          if (flag & I_NO_SEQ) == 0
            total = sum_len + s.l_seq
            mi.ensure_seq_capacity(total)
            s.l_seq.times do |j|
              c = SEQ_NT4_TABLE[s.seq.byte_at(j).to_i].to_i
              mi.seq4_set(sum_len + j, c)
            end
          end
          s.rid = (mi.n_seq - 1).to_u32
          sum_len += s.l_seq
        end

        # Sketch
        seqs.each do |s|
          next if s.l_seq == 0
          Minimap2.mm_sketch(s.seq, s.l_seq, mi.w, mi.k, s.rid, (flag & I_HPC) != 0, p)
        end

        # Add to buckets
        mi.add(p)
        p.clear
      end

      if Minimap2.verbose >= 3
        t = Minimap2.realtime - Minimap2.realtime0
        STDERR.printf("[M::mm_idx_gen::%.3f*%.2f] collected minimizers\n", t, Minimap2.cputime / t)
      end

      mi.post(n_threads)

      if Minimap2.verbose >= 3
        t = Minimap2.realtime - Minimap2.realtime0
        STDERR.printf("[M::mm_idx_gen::%.3f*%.2f] sorted minimizers\n", t, Minimap2.cputime / t)
      end

      mi
    end

    # -------------------------------------------------------------------------
    # Build index from in-memory sequences (mirrors mm_idx_str).
    # -------------------------------------------------------------------------
    def self.from_strings(w : Int32, k : Int32, is_hpc : Bool,
                          bucket_bits : Int32,
                          seqs : Array(String),
                          names : Array(String)? = nil) : MmIdx
      flag = 0
      flag |= I_HPC      if is_hpc
      flag |= I_NO_NAME  if names.nil?
      bucket_bits = 14   if bucket_bits < 0
      mi = new(w, k, bucket_bits, flag)
      mi.n_seq = seqs.size.to_u32
      sum_len = 0_u64
      total_len = seqs.sum(&.size.to_u64)
      if (flag & I_NO_SEQ) == 0
        mi.ensure_seq_capacity(total_len)
      end

      seqs.each_with_index do |s_str, i|
        name = names.try { |ns| ns[i]? } || ""
        offset = sum_len
        len = s_str.size.to_u32
        mi.seq << IdxSeq.new(name, offset, len)
        mi.name_map[name] = i.to_u32 unless name.empty?
        if (flag & I_NO_SEQ) == 0
          len.times do |j|
            c = SEQ_NT4_TABLE[s_str.byte_at(j).to_i].to_i
            mi.seq4_set(offset + j, c)
          end
        end
        if len > 0
          p = [] of Mm128
          Minimap2.mm_sketch(s_str, len.to_i32, w, k, i.to_u32, is_hpc, p)
          mi.add(p)
        end
        sum_len += len
      end

      mi.post(1)
      mi
    end
  end

  # ---------------------------------------------------------------------------
  # Index reader — wraps either a prebuilt .mmi file or a sequence file.
  # Mirrors mm_idx_reader_t / mm_idx_reader_open / mm_idx_reader_read.
  # ---------------------------------------------------------------------------
  class MmIdxReader
    def initialize(@fn : String, @opt : MmIdxOpt, @fn_out : String? = nil)
      @is_idx   = MmIdx.is_idx?(@fn)
      @n_parts  = 0
      if @is_idx > 0
        @idx_io = File.open(@fn, "rb")
        @seq_fp = nil
      else
        @idx_io = nil
        @seq_fp = BSeqFile.open(@fn)
      end
      if fn_out = @fn_out
        @out_io = File.open(fn_out, "wb")
      else
        @out_io = nil
      end
    end

    def read(n_threads : Int32 = 1) : MmIdx?
      mi : MmIdx?
      if (idx = @idx_io)
        mi = MmIdx.load(idx)
        if mi && Minimap2.verbose >= 2
          if mi.k != @opt.k || mi.w != @opt.w || (mi.flag & I_HPC) != (@opt.flag & I_HPC)
            STDERR.puts "[WARNING] Indexing parameters (-k, -w or -H) overridden by parameters used in the prebuilt index."
          end
        end
      elsif fp = @seq_fp
        mi = MmIdx.from_file(
          fp, @opt.w, @opt.k, @opt.bucket_bits, @opt.flag,
          @opt.mini_batch_size.to_i32, n_threads, @opt.batch_size
        )
      end
      if mi
        @out_io.try { |o| mi.dump(o) }
        mi.index = @n_parts
        @n_parts += 1
      end
      mi
    end

    def eof? : Bool
      if @is_idx > 0
        @idx_io.try(&.pos) == @is_idx
      else
        @seq_fp.try(&.eof?) || true
      end
    end

    def close : Nil
      @idx_io.try(&.close)
      @seq_fp.try(&.close)
      @out_io.try(&.close)
    end
  end
end
