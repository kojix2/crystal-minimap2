# Crystal Minimap2

[![build](https://github.com/kojix2/crystal-minimap2/actions/workflows/build.yml/badge.svg)](https://github.com/kojix2/crystal-minimap2/actions/workflows/build.yml)

A crystal implementation of [minimap2](https://github.com/lh3/minimap2).
Created as a practice exercise for Claude Code

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
bin/minimap2 -x map-ont -c ref.fa reads.fa  # with CIGAR
bin/minimap2 -x map-ont -a ref.fa reads.fa  # SAM output
bin/minimap2 -x map-ont -t 8 ref.fa reads.fa
bin/minimap2 -d ref.mmi ref.fa              # build index
bin/minimap2 ref.mmi reads.fa               # map against index
```

Presets: `map-ont`, `map-pb`, `map-hifi`, `asm5`, `asm10`, `asm20`, `splice`, `sr`, `ava-ont`, `ava-pb`

## paftools

```sh
bin/paftools <command> [options] <input>
bin/paftools <command> -h
```

## Limitations

- No paired-end mode
- No split index for very large references
- Single-threaded index loading from `.mmi` files

## Tests

```sh
make spec
```

## License

This is a reimplementation of minimap2. The original license applies to the algorithms and design.
