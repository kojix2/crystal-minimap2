module Paftools
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # paf2gff — convert a protein-to-genome PAF (e.g. miniprot output) to GFF3
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  # Mutable GFF9-ish row: seqname/source/feature/start/end/score/strand/phase/fs.
  # `fs` (frameshift flag, 0/1) is stored in the final attribute column slot
  # until it gets folded into the `frameshift=` attribute string just before
  # printing — mirrors how the JS keeps an Int in a[i][8] until the last loop.
  private class GRow
    property seqname : String
    property source : String
    property feature : String
    property st : Int32
    property en : Int32
    property score : String
    property strand : String
    property phase : String
    property fs : Int32

    def initialize(@seqname, @source, @feature, @st, @en, @score, @strand, @phase, @fs)
    end
  end

  def self.cmd_paf2gff(args : Array(String)) : Int32
    aa = false
    rest = [] of String
    i = 0
    while i < args.size
      case args[i]
      when "-a"          ; aa = true
      when "-h", "--help"; STDERR.puts "Usage: paftools paf2gff [-a] <in.paf>"; return 0
      else                 rest << args[i]
      end
      i += 1
    end
    if rest.empty?
      STDERR.puts "Usage: paftools paf2gff [-a] <in.paf>"
      return 1
    end

    re_cigar = /(\d+)([A-Z=])/
    hid = 1; last_name : String? = nil

    open_in(rest[0]) do |io|
      io.each_line(chomp: true) do |line|
        t = line.split('\t')
        next if t[5] == "*" # skip unmapped lines

        if t[0] != last_name
          last_name = t[0]; hid = 1
        else
          hid += 1
        end

        qlen = t[1].to_i; qs = t[2].to_i; qe = t[3].to_i
        tlen = t[6].to_i; ts = t[7].to_i; te = t[8].to_i
        nmatch = t[9].to_i; blen = t[10].to_i
        strand = t[4]

        cigar : String? = nil; score : Int32? = nil; np : Int32? = nil
        dist_stop : Int32? = nil; dist_start : Int32? = nil
        (12...t.size).each do |j|
          if m = /^(cg:Z|AS:i|np:i|da:i|do:i):(\S+)/.match(t[j])
            case m[1]
            when "cg:Z"; cigar = m[2]
            when "AS:i"; score = m[2].to_i
            when "np:i"; np = m[2].to_i
            when "do:i"; dist_stop = m[2].to_i
            when "da:i"; dist_start = m[2].to_i
            end
          end
        end
        raise "failed to find the cg:Z tag" unless cigar
        raise "failed to find the AS:i tag" unless score

        st = 0; en = 0; phase = 0; pseudo = false; fs = 0
        a = [] of GRow
        if dist_start && dist_start == 0
          a << GRow.new(t[5], "paf2gff", "start_codon", 0, 3, "0", strand, ".", 0)
        end
        cigar.scan(re_cigar) do |m|
          len = m[1].to_i
          op = m[2]
          case op
          when "M", "D"
            en += aa ? len * 3 : len
          when "F", "G", "R"
            en += len; pseudo = true; fs = 1
          when "N"
            a << GRow.new(t[5], "paf2gff", "exon", st, en, "0", strand, phase.to_s, fs)
            st = en + len; en += len; phase = 0; fs = 0
          when "U" # ...xGT...AGxx...
            a << GRow.new(t[5], "paf2gff", "exon", st, en + 1, "0", strand, phase.to_s, fs)
            st = en + len - 2; en += len; phase = 2; fs = 0
          when "V" # ...xxGT...AGx...
            a << GRow.new(t[5], "paf2gff", "exon", st, en + 2, "0", strand, phase.to_s, fs)
            st = en + len - 1; en += len; phase = 1; fs = 0
          end
        end
        a << GRow.new(t[5], "paf2gff", "exon", st, en, "0", strand, phase.to_s, fs)
        raise "inconsistent cigar" if en != te - ts
        if dist_stop && dist_stop == 0
          a << GRow.new(t[5], "paf2gff", "stop_codon", en, en + 3, "0", strand, ".", 0)
        end

        type = pseudo ? "pseudogene" : "protein_coding"
        attr = "transcript_id=#{t[0]}##{hid};transcript_type=#{type}"
        # js_fixed matches JS toFixed's "NaN"/"Infinity" output for the
        # (unlikely but possible) blen==0 edge case, instead of Crystal
        # sprintf's "nan"/"inf".
        trans_attr = "identity=#{js_fixed(nmatch.to_f / blen, 4)}"
        if npv = np
          trans_attr += ";positive=#{js_fixed(npv * 3.0 / blen, 4)}"
        end
        trans_attr += ";aa_start=#{qs}"
        trans_attr += ";aa_end=#{qlen - qe}"
        trans_attr += ";dist_start_codon=#{dist_start}" if dist_start && dist_start >= 0
        trans_attr += ";dist_stop_codon=#{dist_stop}" if dist_stop && dist_stop >= 0
        trans_st = ts; trans_en = te
        if dist_stop && dist_stop == 0
          if strand == "-"
            trans_st -= 3
          else
            trans_en += 3
          end
        end
        puts [t[5], "paf2gff", "transcript", trans_st + 1, trans_en, score, strand, ".",
              attr + ";" + trans_attr].join('\t')

        if aa && strand == "-"
          len = te - ts
          b = [] of GRow
          (a.size - 1).downto(0) do |k|
            row = a[k]
            x = len - row.st
            row.st = len - row.en
            row.en = x
            # JS keeps this line commented out — "not sure if this line is
            # needed" — original author's own uncertainty, so we leave the
            # phase field untouched here too rather than guess:
            # row.phase = row.phase == "0" ? "0" : (3 - row.phase.to_i).to_s
            b << row
          end
          a = b
        end

        a.each do |row|
          row.feature = "CDS" if !pseudo && row.feature == "exon"
          row.st += ts + 1
          row.en += ts
          puts [row.seqname, row.source, row.feature, row.st, row.en, row.score,
                row.strand, row.phase, "#{attr};frameshift=#{row.fs}"].join('\t')
        end
      end
    end
    0
  end
end
