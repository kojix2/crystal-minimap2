module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # pafcmp — compare a "test" PAF against a "base" (truth) PAF, per query name
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  private struct PafcmpBase
    property tname : String
    property ts : Int32
    property te : Int32
    property mapq : Int32
    property n_hit : Int32
    property n_wrong : Int32

    def initialize(@tname, @ts, @te, @mapq, @n_hit = 0, @n_wrong = 0)
    end
  end

  def self.cmd_pafcmp(args : Array(String)) : Int32
    min_len = 5000; min_mapq_opt = 10; min_ovlp = 0.5
    rest = [] of String; i = 0
    while i < args.size
      case args[i]
      when "-q"          ; i += 1; min_mapq_opt = args[i].to_i
      when "-h", "--help"; STDERR.puts "Usage: paftools pafcmp [options] <base.paf> <test.paf>"; return 0
      else                 rest << args[i]
      end
      i += 1
    end
    if rest.size < 2
      STDERR.puts "Usage: paftools pafcmp [options] <base.paf> <test.paf>"
      STDERR.puts "Options:"
      STDERR.puts "  -q INT    min mapping quality [#{min_mapq_opt}]"
      return 1
    end

    n_base = 0; n_test = 0
    # NOTE (bug fix): the JS original increments `opt.n_out_high`/
    # `opt.n_out_low` here — but those counters were declared on the `eval`
    # object (`n_out_high:0, n_out_low:0`), not on `opt`. Because of that typo
    # the real eval.n_out_high/n_out_low are *always* 0, so the final
    # "additional test alignments" line always prints 0 regardless of actual
    # unmatched test alignments. We track these correctly here instead of
    # replicating the no-op counters.
    n_out_high = 0; n_out_low = 0
    n_hit = 0; n_wrong = 0; n_miss = 0

    base = Hash(String, PafcmpBase).new

    process_base = ->(a : Array(Array(String))) {
      if a.size == 1
        row = a[0]
        blen = row[1].to_i
        if blen >= min_len
          mapq = row[11].to_i
          n_base += 1 if mapq >= min_mapq_opt
          base[row[0]] = PafcmpBase.new(row[5], row[7].to_i, row[8].to_i, mapq)
        end
      end
    }

    open_in(rest[0]) do |io|
      STDERR.puts "Reading #{rest[0]}..."
      a = [] of Array(String)
      io.each_line(chomp: true) do |line|
        next if /\ttp:A:S/.matches?(line)
        t = line.split('\t')
        if a.size > 0 && a[0][0] != t[0]
          process_base.call(a)
          a = [] of Array(String)
        end
        a << t
      end
      process_base.call(a)
    end

    process_test = ->(a : Array(Array(String))) {
      row = a[0]
      blen = row[1].to_i
      if blen >= min_len
      mapq = row[11].to_i
      n_test += 1 if mapq >= min_mapq_opt
      c_tname = row[5]; c_ts = row[7].to_i; c_te = row[8].to_i; c_mapq = mapq
      b = base[row[0]]?
      if b.nil?
        if c_mapq >= min_mapq_opt
          n_out_high += 1
        else
          n_out_low += 1
        end
      else
        inter = 0
        union = (b.te - b.ts) + (c_te - c_ts)
        if b.tname == c_tname # same chr
          if b.ts < c_ts
            if b.te > c_ts
              inter = b.te - c_ts
              union = c_te - b.ts
            end
          else # c.ts < b.ts
            if c_te > b.ts
              inter = c_te - b.ts
              union = b.te - c_ts
            end
          end
        end
        if inter >= union * min_ovlp
          n_hit += 1 if b.mapq >= min_mapq_opt
          b.n_hit += 1
        else
          if b.mapq >= min_mapq_opt
            puts ["W", row[0], b.tname, b.ts, b.te, b.mapq, c_tname, c_ts, c_te, c_mapq].join('\t')
            n_wrong += 1
          end
          b.n_wrong += 1
        end
        base[row[0]] = b
      end
      end
    }

    open_in(rest[1]) do |io|
      STDERR.puts "Reading #{rest[1]}..."
      a = [] of Array(String)
      io.each_line(chomp: true) do |line|
        next if /\ttp:A:S/.matches?(line)
        t = line.split('\t')
        if a.size > 0 && a[0][0] != t[0]
          process_test.call(a)
          a = [] of Array(String)
        end
        a << t
      end
      process_test.call(a) unless a.empty?
    end

    base.each do |r, b|
      if b.mapq >= min_mapq_opt && b.n_hit == 0 && b.n_wrong == 0
        n_miss += 1
        puts ["M", r, b.tname, b.ts, b.te, b.mapq].join('\t')
      end
    end

    puts "X\t#{n_base} base alignments with mapQ>=#{min_mapq_opt}"
    puts "X\t#{n_hit} base alignments correctly mapped by test"
    puts "X\t#{n_wrong} wrong test alignment"
    puts "X\t#{n_miss} base alignments missing"
    puts "X\t#{n_out_high} additional test alignments with mapQ>=#{min_mapq_opt}"
    0
  end
end
