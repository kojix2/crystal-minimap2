# Crystal Minimap2

[![build](https://github.com/kojix2/crystal-minimap2/actions/workflows/build.yml/badge.svg)](https://github.com/kojix2/crystal-minimap2/actions/workflows/build.yml)

A complete Crystal implementation of [minimap2](https://github.com/lh3/minimap2) — a versatile pairwise sequence aligner for nucleotide sequences. All algorithms are implemented directly in Crystal with no C bindings.

## Requirements

- Crystal >= 1.19.1

## Build

```sh
make
```

Builds both `bin/minimap2` and `bin/paftools`. Individual targets:

```sh
make minimap2
make paftools
```

## Usage

```sh
bin/minimap2 -x map-ont ref.fa reads.fa
bin/minimap2 -x map-ont -c ref.fa reads.fa     # with CIGAR in PAF
bin/minimap2 -x map-ont -a ref.fa reads.fa     # SAM output
bin/minimap2 -x map-ont --cs ref.fa reads.fa   # cs tag (short form)
bin/minimap2 -x map-ont --MD -a ref.fa reads.fa # SAM with MD tag
bin/minimap2 -x map-ont -t 8 ref.fa reads.fa  # 8 threads
bin/minimap2 -d ref.mmi ref.fa                 # build index
bin/minimap2 ref.mmi reads.fa                  # map against index
```

Presets: `map-ont`, `map-pb`, `map-hifi`, `asm5`, `asm10`, `asm20`, `splice`, `splice:hq`, `splice:sr`, `sr`, `ava-ont`, `ava-pb`

### Key options

| Option | Description |
|--------|-------------|
| `-x STR` | Preset (see above) |
| `-a` | SAM output |
| `-c` | CIGAR in PAF |
| `--cs` | cs tag (short form) |
| `--cs-long` | cs tag (long form) |
| `--MD` | MD tag |
| `--eqx` | =/X CIGAR operators |
| `-t INT` | Number of threads |
| `-N INT` | Max secondary alignments |
| `-G NUM` | Max intron length |
| `-C INT` | Non-canonical splice penalty |
| `-u CHAR` | GT-AG direction (f/b/n/r) |
| `--secondary=no` | Suppress secondary alignments |

## paftools

```sh
bin/paftools <command> [options] <input>
bin/paftools <command> -h
```

## Limitations

- No SIMD acceleration (scalar int32 alignment; slower than C for large alignments)
- No split index for very large references (index is built in memory)
- Single-threaded index loading from `.mmi` files

## Tests

```sh
make spec
```

## License

This is a reimplementation of minimap2. The original license applies to the algorithms and design.
