# crystal-minimap2

A pure-Crystal port of [minimap2](https://github.com/lh3/minimap2) — no C bindings, no FFI.
All core algorithms are re-implemented in Crystal.

## Requirements

- Crystal >= 1.19.1

## Build

```sh
shards build --release
```

The binary is written to `bin/minimap2`.
Multi-threading uses `Fiber::ExecutionContext::Parallel` and is enabled automatically by the build flags in `shard.yml`.

## Usage

```sh
# Map long reads (PAF output)
bin/minimap2 -x map-ont ref.fa reads.fa

# Map with CIGAR
bin/minimap2 -x map-ont -c ref.fa reads.fa

# SAM output
bin/minimap2 -x map-ont -a ref.fa reads.fa

# Use multiple threads
bin/minimap2 -x map-ont -t 8 ref.fa reads.fa

# Build a prebuilt index, then map
bin/minimap2 -d ref.mmi ref.fa
bin/minimap2 ref.mmi reads.fa

# Splice-aware alignment
bin/minimap2 -x splice genome.fa rna.fa
```

## Presets

`map-ont`, `map-pb`, `map-hifi`, `asm5`, `asm10`, `asm20`, `splice`, `sr`, `ava-ont`, `ava-pb`

## Limitations

Compared to the C minimap2:

- No paired-end mode (`pe.c`)
- No split index for very large references (`splitidx.c`)
- No `paftools.js` utilities
- Single-threaded index loading from `.mmi` files (mapping itself is parallel)

## Tests

```sh
crystal spec -Dpreview_mt -Dexecution_context
```

## License

MIT
