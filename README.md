# hifi_trimmer

`hifi-trimmer` is a command-line tool for filtering and trimming extraneous adapter hits
from a HiFi read set using a BLAST search against a fasta file of adapter sequences. It is
designed to be highly configurable, with per-adapter settings to determine actions if
the adapter is found at the ends of a read or in the middle. To improve reproducibility,
the primary output of the tool is a BED file that describes the region of each read to be
excluded. The tool also includes a command to filter the reads to disk using the produced
BED file.

The polars backend for BLAST file processing should respect the number of cores set
by your scheduler; however, if this is not the case the number of threads used can be
adjusted by setting the environment variable `POLARS_MAX_THREADS=int` before running the
software. The number of threads used by the  bgzip backend for writing compressed
BED files and the filtered FASTA files can be adjusted using the command-line option
`--threads`.

## Installation

```
pip install hifi_trimmer
```

## Usage

```
Usage: hifi-trimmer [OPTIONS] COMMAND [ARGS]...

  Main entry point for the tool.

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  process-blast  Processes the input blastout file according to the...
  trim           Filter the reads stored in a BAM file using the...
```

To process a blast output TSV to a BED file:

```
Usage: hifi-trimmer process-blast [OPTIONS] BLASTOUT ADAPTER_YAML

  Processes the input blastout file according to the adapter yaml key.

  BLASTOUT: tabular file resulting from a BLAST query of a readset against a
  BLAST database of adapter sequences, run with -outfmt "6 std qlen".

  ADAPTER_YAML: yaml file containing a list with the following fields per
  adapter:

  - name: (name of adapter. can be a regular expression)
    discard_middle: True/False (discard read if adapter found in middle)
    discard_end: True/False (discard read if adapter found in end)
    trim_end: True/False (trim read if adapter found in end)
    middle_pident: int (minimum pident requred to identify adapter in middle of read)
    middle_length: int (minimum match length required to identify adapter in middle of read)
    end_pident: int (minimum pident requred to identify adapter in end window)
    end_length: int (minimum match length requred to identify adapter in end window)

  Output: By default, writes bgzipped BED to [prefix].bed.gz, and a JSON
  summary file with raw counts of adapter hits detected, counts identified
  after processing, and the total length of removed sequences per adapter to
  [prefix].summary.json.

Options:
  -p, --prefix TEXT               Output prefix for results. Defaults to the
                                  basename of the blastout if not provided.
  -ml, --min-length-after-trimming INTEGER
                                  Minumum length of a read after trimming the
                                  ends in order not to be discarded  [default:
                                  300]
  -el, --end-length INTEGER       Window size at either end of the read to be
                                  considered as 'ends' for searching
                                  [default: 150]
  --hits / --no-hits              Write the hits identified using the given
                                  adapter specifications to TSV. The format is
                                  standard BLAST outfmt 6 with the following
                                  extra columns: read_length (int), discard
                                  (bool), trim_l (bool), trim_r (bool)
  --no-summary                    Skip writing a summary JSON with the number
                                  of hits for each adapter
  -t, --threads INTEGER           Number of threads to use for compression
                                  [default: 1]
  --help                          Show this message and exit.
```

To filter a BAM file using the BED file:

```
Usage: hifi-trimmer trim [OPTIONS] BAM BED

  Filter the reads stored in a BAM file using the appropriate BED file
  produced by blastout_to_bed and write to a bgzipped fasta file.

  BAM: BAM file in which to filter reads
  BED: BED file describing regions of the read set to exclude.
  OUTFILE: Output trimmed BAM file. Defaults to stdout.

Options:
  -o, --outfile FILENAME
  -t, --threads INTEGER        Number of threads to use for compression
                               [default: 1]
  -f, --format [sam|bam|cram]  Output file format.  [default: bam]
  -h, --format-opt TEXT        Format options to pass to htslib for writing
                               the output (see https://www.htslib.org/doc/samt
                               ools.html#GLOBAL_COMMAND_OPTIONS). Can be
                               specified multiple times if multiple options
                               are required.
  --help                       Show this message and exit.
```

## Example

First, BLAST your reads against an adapter database:

```
blastn -query <(samtools fasta /path/to/bam) \
  -db /path/to/adapter/blast/db \
  -reward 1 -penalty -5 -gapopen 3 -gapextend 3 \
  -dust no -soft_masking true -evalue 700 \
  -searchsp 1750000000000 \
  -outfmt "6 std qlen" |\
  bgzip > blastout.gz
```

To create a BED file, you then need to create a YAML file describing the actions to take for each adapter. The adapter name
can be a regular expression, but note that each adapter name in the BLAST file must match only one entry in
the YAML.

```
- adapter: "^NGB00972"  // regular expression matching adapter names
  discard_middle: True  // discard read if adapter is found in the middle
  discard_end: False    // discard read if adapter found at end
  trim_end: True        // trim read if adapter is found at end (overridden by discard choice)
  middle_pident: 95     // minimum percent identity for a match in the middle of the read
  middle_length: 44     // minimum match length for a match in the middle of the read
  end_pident: 90        // minimum percent identity for a match at the end of the read
  end_length: 18        // minimum match length for a match at the end of the read
```

Then run `process-blast` to generate a BED file:

```
hifi-trimmer process-blast /path/to/blastout.gz /path/to/yaml
```

Then filter the BAM file using the BED file:

```
hifi-trimmer trim /path/to/bam /path/to/bed /path/to/final/file.bam
```
