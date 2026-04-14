require "./minimap2/constants"
require "./minimap2/types"
require "./minimap2/misc"
require "./minimap2/sketch"
require "./minimap2/sdust"
require "./minimap2/options"
require "./minimap2/bseq"
require "./minimap2/index"
require "./minimap2/seed"
require "./minimap2/lchain"
require "./minimap2/ksw2"
require "./minimap2/hit"
require "./minimap2/align"
require "./minimap2/format"
require "./minimap2/map"
require "./minimap2/cli"

# Pure-Crystal port of [minimap2](https://github.com/lh3/minimap2),
# a versatile pairwise sequence aligner for nucleotide sequences.
#
# All algorithms are implemented directly in Crystal — no C bindings.
#
# Quick start:
# ```
# aligner = Minimap2::Aligner.from_strings(["ACGT" * 100], ["ref"], "map-ont")
# aligner.map("ACGTACGT...") do |hit|
#   puts "#{hit.rid}\t#{hit.rs}-#{hit.re}\t#{hit.rev ? '-' : '+'}"
# end
# ```
module Minimap2
  VERSION = "0.1.0"
end
