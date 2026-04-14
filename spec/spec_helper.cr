require "spec"
require "../src/minimap2"

# Suppress minimap2 progress messages during tests
# Keep only error output during tests (suppress [M::] progress messages)
Minimap2.verbose = 2
