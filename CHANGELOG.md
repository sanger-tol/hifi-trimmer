# hifi_trimmer: changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[1.2.2](https://github.com/sanger-tol/hifi-trimmer/releases/tag/v1.2.2)] - [2025-07-09]

### Fixed
- Remove extraneous duckdb dependency

## [[1.2.1](https://github.com/sanger-tol/hifi-trimmer/releases/tag/v1.2.1)] - [2025-07-09]

### Fixed
- Fixes a bug where reads which should be discarded due to an end adapter are skipped.

## [[1.2.0](https://github.com/sanger-tol/hifi-trimmer/releases/tag/v1.2.0)] - [2025-07-09]

### Added
- Adds a `--fastq` option to the `filter_bam` command, which writes the trimmed reads in FASTQ format instead of FASTA

## [[1.1.0](https://github.com/sanger-tol/hifi-trimmer/releases/tag/v1.1.0)] - [2025-02-28]

### Added
- `--version` flag
- Default values are included in the `--help` output for options.

### Fixed
- Rare problem encountered where read lengths were over 65,535 - internally, read length and related variables are now typed as UInt32 rather than UInt16.
- Refactored code in the `determine_actions` function to be more readable and to avoid intermittent errors with typing.

## [[1.0.1](https://github.com/sanger-tol/hifi-trimmer/releases/tag/v1.0.1)] - [2025-02-28]

### Fixed
- Fix rare bug where for some reason the Polars StringCache ends up non-comparable by enabling the global StringCache.
- This reduces performance slightly but it should not be noticable in most cases, unless there are hundreds of millions of rows in the BLAST file.

## [[1.0.0](https://github.com/sanger-tol/hifi-trimmer/releases/tag/v1.0.0)] - [2025-02-27]

Base release of hifi_trimmer.

### Added
- process_blast: Process a BLAST output file and determine which reads to discard or trim, producing a compliant BED file
- filter_bam: Process a BAM file using the BED file to produce a filtered FASTA file
