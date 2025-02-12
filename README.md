# hifi_trimmer

## Installation

```
git clone git@github.com:prototaxites/hifi-trimmer.git
cd hifi-trimmer
pip install .
```

## Usage

```
$ uv run hifi_trimmer blastout_to_bed --help
Usage: hifi_trimmer blastout_to_bed [OPTIONS] BLASTOUT ADAPTER_YAML

  Processes the input blastout file according to the adapter yaml key.

  BLASTOUT: tabular file resulting from BLAST with -outfmt "6 std qlen". If
  the qlen column is missing, lengths can be calculated by passing the --bam
  option. ADAPTER_YAML: yaml file contaning a list with the following fields
  per adapters:     - name: name of adapter. can be a regular expression
  discard_middle: True/False - discard read if adapter found in middle
  discard_end: True/False  - discard read if adapter found in end
  trim_end: True/False - trim read if adapter found in end
  middle_pident: int - minimum pident requred to identify adapter in middle of
  read       middle_length: int - minimum match length required to identify
  adapter in middle of read       end_pident: int - minimum pident requred to
  identify adapter in end window       end_length: int - minimum match length
  requred to identify adapter in end window

  Output: BED file to stdout

Options:
  --bam PATH                      If blastout file has no read length field, a
                                  BAM file of reads to get read lengths
  --min_length_after_trimming INTEGER
                                  Minumum length of a read after trimming the
                                  ends in order not to be discarded.
  --end_length INTEGER            Window size at either end of the read to be
                                  considered as 'ends' for searching.
  --help                          Show this message and exit.
```

```
$ uv run hifi_trimmer filter_bam_to_fasta --help
Usage: hifi_trimmer filter_bam_to_fasta [OPTIONS] BED BAM [OUTFILE]

  Filter the reads stored in a BAM file using the appropriate BED file
  produced  by blastout_to_bed and write to a bgzipped fasta file.

  BED: BED file describing regions of the read set to exclude. BAM: BAM file
  in which to filter reads OUTFILE: (optional) The output bgzipped fasta file

Options:
  --help  Show this message and exit.
```

## Example

First, BLAST your reads against an adapter database:

```
blast -query <(samtools fasta /path/to/bam) -db /path/to/adapter/blast/db -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust no -soft_masking true -evalue 700 -searchsp 1750000000000 -outfmt "6 std qlen" | bgzip > blastout.gz
```

Then, create a YAML file describing what to do with each adapter hit:
```
- adapter: "NGB00972"
  discard_middle: True // Discard read if hit found in middle of read
  discard_end: False // Discard read if hit found in end of read
  trim_end: True // Trim read if hit found at end of read and read not otherwise discarded
  middle_pident: 95 // Minimum hit percent identity for match in middle
  middle_length: 44 // Minimum hit length for match in middle
  end_pident: 90 // Minimum hit percent identity for match at end
  end_length: 18  // Minimum hit length for match at end
- adapter: "NGB00973"
  discard_middle: True
  discard_end: False
  trim_end: True
  middle_pident: 95
  middle_length: 34
  end_pident: 90
  end_length: 18
```

Using these, you can process the BLAST results into a BED file:

```
hifi_trimmer blastout_to_bed blastout.gz adapters.yaml > out.bed
```

And then process the BAM and BED files:

```
hifi_trimmer filter_bam_to_fasta out.bed /path/to/bam reads.filtered.fa.gz
```