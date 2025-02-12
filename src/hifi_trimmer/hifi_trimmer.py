import click
import bgzip
from pathlib import Path
import polars as pl
import polars.selectors as cs
import pysam
import re
import sys
import csv
import yaml

def read_adapter_yaml(yaml_path: str) -> pl.DataFrame:
    """Open an adapter YAML file and return it as a lazy pl.DataFrame"""
    with open(yaml_path) as stream:
        try:
            adapters = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    return pl.DataFrame(adapters).lazy()

def read_blast(blast_path: str) -> pl.LazyFrame:
    ## Read BLAST as lazy DataFrame - will only be streamed from file after 
    ## .collect() called
    blastout = pl.scan_csv(blast_path, has_header=False, separator='\t')

    ## rename blast columns based on number of input columns
    blast_ncol = len(blastout.head(1).collect().columns)

    column_names = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]
    if blast_ncol == 13:
        column_names = column_names + ["read_length"]

    blastout = (blastout
        .rename({old: new for old, new in zip(blastout.collect_schema().names(), column_names)})
        ## reorder start-end to correct strandedness as we're not interested here
        .with_columns(
            pl.when(
                pl.col("qstart") > pl.col("qend")
            )
            .then(pl.struct(
                start = "qend",
                end = "qstart",
            ))
            .otherwise(pl.struct(["qstart", "qend"]))
            .struct.field(["qstart", "qend"])
        )
    )
    
    return blastout

def read_lengths(filtered_reads: set, bam: str) -> pl.DataFrame:
    """Count the lengths of reads in set filtered_reads using bamfile bam"""
    lengths = []
    with pysam.AlignmentFile(bam, "rb", check_sq=False, require_index=False) as b:
        for read in b.fetch(until_eof=True):
            if read.query_name in filtered_reads:
                lengths.append({
                    'qseqid': read.query_name,
                    'read_length': read.query_length
                })

    return pl.DataFrame(lengths)

def check_records(blast: pl.LazyFrame, adapters: pl.DataFrame) -> bool:
    """Check if any adapters in a BLAST dataframe match to more than one
    adapter sequence in the YAML file.    
    """
    adapters = (adapters.select("adapter").unique())
    counts = (blast
        .select(pl.col("sseqid"))
        .unique()
        .join(adapters, how = "cross")
        .filter(
            pl.struct(["sseqid", "adapter"]).map_elements(
                lambda s: bool(re.search(s["adapter"], s["sseqid"])), return_dtype=pl.Boolean
            )
        )
        .group_by("sseqid")
        .agg(n = pl.col("sseqid").n_unique())
    ).collect()["n"].to_list()

    return any(x > 1 for x in counts)

def format_fasta_record(header: str, sequence: str) -> str:
    """Format a header and sequence into FASTA format"""
    return f">{header}\n{sequence}\n"

def iter_exhausted(iter):
    """Check if the provided iterable is exhausted."""
    try:
        next(iter)
        return False
    except StopIteration:
        return True

def trim_positions(seq, ranges):
    """Trim DNA sequence seq to remove the positions specified in ranges.
    
    seq: string
    ranges: list of non-overlapping tuples with (start, end)
    """
    ## sort the ranges so we trim from the end - that way indexing stays constant
    ranges = sorted(ranges, key=lambda x: x[0], reverse=True)
    
    for start, end in ranges:
        seq = seq[:start] + seq[end:]
    
    return seq

def match_hits(blastout: pl.DataFrame, adapter_df: str, end_length: int) -> pl.DataFrame:
    """Determine whether each blast hit matches sufficient criteria for discard/trimming
    
    Keyword arguments:
    blastout: pl.DataFrame containing the results of a blast search with read lengths in column 13
    adapter_key: pl.DataFrame of adapters and their filtering/trimming parameters
    end_length: length of window at either end of read which is examined for trimming
    """
    ## When to discard a read due to adapter presence
    discard = (
        ((pl.col("discard_middle")).and_(
            (pl.col("read_length") > 2 * end_length),
            (pl.col("qstart") > end_length),
            (pl.col("qend") < pl.col("read_length") - end_length),
            (pl.col("pident") >= pl.col("middle_pident")),
            (pl.col("length") >= pl.col("middle_length"))
        )) |
        ((pl.col("discard_end")).and_(
            (pl.col("qend") < end_length),
            (pl.col("qstart") > pl.col("read_length") - end_length),
            (pl.col("pident") >= pl.col("end_pident")),
            (pl.col("length") >= pl.col("end_length"))
        ))
    )

    ## When to trim a read at the left end
    trim_l = (
        (pl.col("trim_end"))
        .and_(
            pl.col("qend") <= end_length,
            (pl.col("pident") >= pl.col("end_pident")),
            (pl.col("length") >= pl.col("end_length"))
        )
    )

    ## When to trim a read at the righr end
    trim_r = (
        (pl.col("trim_end"))
        .and_(
            pl.col("qstart") >= pl.col("read_length") - end_length,
            (pl.col("pident") >= pl.col("end_pident")),
            (pl.col("length") >= pl.col("end_length"))
        )
    )
    
    blastout = (blastout
        ## match all adapters against all rows and filter to keep just those that
        ## have a regex match
        .join(adapter_df, how="cross", on=None, maintain_order="left")
        .filter(
            pl.struct(["sseqid", "adapter"]).map_elements(
                lambda s: bool(re.search(s["adapter"], s["sseqid"])), return_dtype=pl.Boolean
            )
        )
        ## Determine which hits to keep and return rows with a "real" hit
        .with_columns(
            discard=discard,
            trim_l=trim_l,
            trim_r=trim_r
        )
        .filter(pl.any_horizontal("discard", "trim_l", "trim_r"))
    )

    return blastout

def create_bed(blastout: pl.DataFrame, end_length: int, min_length: int) -> pl.DataFrame:
    """Take a dataframe with determined matches from match_hits()
    and process the result into a bedfile.
    
    Keyword arguments:
    blastout: pl.DataFrame of blast hits with matches
    end_length: length of window at either end of read which is examined for trimming
    min_length: minimum length of a read after trimming to keep
    """
    ## Calculate read length after trimming (naiively) and then unpivot the trim columns so the
    ## trim section below operates separately for l and r
    blastout = (blastout
        .with_columns(
            (pl.col("read_length") - ((pl.col("trim_l").any() + pl.col("trim_r").any()) * end_length)).over("qseqid").alias("read_length_after_trimming"),
        )
        .unpivot(
            cs.matches("trim_(l|r)"),
            index=~cs.matches("trim_(l|r)"),
            variable_name="trim_where",
            value_name="trim"
        )
    )

    ## Identify reads to discard and format as BED file with one line
    ## for the whole read
    discard = (blastout
        .filter((pl.col("discard").any() | pl.col("read_length_after_trimming").lt(min_length).any()).over("qseqid"))
        .group_by("qseqid")
        .agg(
            start = pl.lit(0, dtype=pl.Int64),
            end = pl.col("read_length").first(),
            reason = pl.concat_str(pl.lit("discard:"), pl.col("sseqid").first())
        )
    )

    ## Identify reads that are not to be discarded and then build a BED file
    ## for left, right, or both depending on what is present
    trim = (
        blastout
        ## Remove discarded reads
        .filter(~(pl.col("discard").any() | pl.col("read_length_after_trimming").lt(min_length).any()).over("qseqid"))
        ## Filter so that we only keep rows where we are trimming
        .filter(pl.col("trim"))
        ## Using structs, return start and end trim positions depending on whether
        ## left trim or right - and convert to columns
        .with_columns(
            pl.when(
                (pl.col("trim_where") == "trim_l")
            ).then(pl.struct(
                pl.lit(0).alias("start"),
                pl.lit(end_length).alias("end")
                )
            ).otherwise(pl.struct(
                (pl.col("read_length") - end_length).alias("start"),
                pl.col("read_length").alias("end")
                )
            ).struct.field(["start", "end"])
        )
        ## Format as BED
        .group_by("qseqid", "trim_where")
        .agg(
            pl.col("start").min().alias("start"),
            pl.col("end").max().alias("end"),
            pl.concat_str(
                pl.lit("trim:"),
                pl.col("sseqid").get(pl.col("qend").arg_max())
            ).alias("reason")
        )
        .select(["qseqid", "start", "end", "reason"])
    )

    ## Concatenate discard and trim and return
    result = pl.concat([discard, trim])
    return result

@click.group()
def cli():
    """Main entry point for the tool."""
    pass

@click.command("blastout_to_bed")
@click.argument('blastout', type=click.Path(exists=True))
@click.argument('adapter_yaml', type=click.Path(exists=True))
@click.option('--bam', default=None, required=False, type=click.Path(exists=True), 
    help="If blastout file has no read length field, a BAM file of reads to get read lengths")
@click.option("--min_length_after_trimming", default=300, type=int,
    help="Minumum length of a read after trimming the ends in order not to be discarded.")
@click.option("--end_length", default=150, type=int,
    help = "Window size at either end of the read to be considered as 'ends' for searching.")
def blastout_to_bed(blastout, adapter_yaml, bam, min_length_after_trimming, end_length):
    """Processes the input blastout file according to the adapter yaml key.
    
    BLASTOUT: tabular file resulting from BLAST with -outfmt "6 std qlen". If the qlen column
    is missing, lengths can be calculated by passing the --bam option.
    ADAPTER_YAML: yaml file contaning a list with the following fields per adapters:
        - name: name of adapter. can be a regular expression
          discard_middle: True/False - discard read if adapter found in middle
          discard_end: True/False  - discard read if adapter found in end
          trim_end: True/False - trim read if adapter found in end
          middle_pident: int - minimum pident requred to identify adapter in middle of read
          middle_length: int - minimum match length required to identify adapter in middle of read
          end_pident: int - minimum pident requred to identify adapter in end window
          end_length: int - minimum match length requred to identify adapter in end window

    Output:
    BED file to stdout
    """ 
    adapters = read_adapter_yaml(adapter_yaml)
    blast = read_blast(blastout)

    ## Make sure each barcode that matches an adapter matches only one adapter
    if check_records(blast, adapters):
        sys.exit(f"Error: adapters in {adapter_yaml} match to more than one adapter in {blastout}!")
    
    ## if read lengths not provided but bam is - get them from the bam
    if(len(blast.head(1).collect().columns) == 12 and bam is not None):
        filtered_reads = set(blast.select("qseqid").unique().collect()["qseqid"].to_list())
        lengths = read_lengths(filtered_reads, bam)
        blast = blast.join(lengths.lazy(), on = "qseqid", how = "left", suffix="")
    else:
        sys.exit("Error: BLAST file only has 12 columns, but no BAM file has been provided to count!")
    
    ## process the blastout file
    hits = match_hits(blast, adapters, end_length)
    bed = (create_bed(hits, end_length, min_length_after_trimming)
        .sort(by=[pl.col("qseqid").str.extract(r"/(\d+)/").cast(pl.Int64), "start"])
        .collect()
    )

    ## Print to stdout
    click.echo(bed.write_csv(separator="\t", include_header=False))

@click.command("filter_bam_to_fasta")
@click.argument('bed', type=click.Path(exists=True))
@click.argument('bam', type=click.Path(exists=True))
@click.argument('outfile', default=None, required=False, type=click.Path(exists=False))
def filter_bam_to_fasta(bed, bam, outfile):
    """
    Filter the reads stored in a BAM file using the appropriate BED file produced 
    by blastout_to_bed and write to a bgzipped fasta file.

    BED: BED file describing regions of the read set to exclude.
    BAM: BAM file in which to filter reads
    OUTFILE: (optional) The output bgzipped fasta file
    """
    if outfile is None:
        outfile = Path(bam).stem + ".filtered.fa.gz"

    ## read BED as dictionary
    filters = csv.DictReader(open(bed, "r"), delimiter="\t", fieldnames=["read", "start", "end", "reason"])
    with open(outfile, "wb") as raw:
        with bgzip.BGZipWriter(raw) as out:
            with pysam.AlignmentFile(bam, "rb", check_sq=False, require_index=False) as b:
                
                ## initialise first BED record
                r = next(filters, None)
                ## Process: for each read in the BAM, check if it matches the current BED record.
                ## If yes, pull BED records until we reach a record for the next read.
                ## Then trim the sequence using the records pulled and write to fasta.
                ## If not, write record straight to disk
                for read in b.fetch(until_eof=True):
                    if r is not None and r["read"] == read.query_name:
                        ranges = [(int(r["start"]), int(r["end"]))]

                        while True:
                            r = next(filters, None)
                            if r is None or r["read"] != read.query_name:
                                break
                            ranges.append((int(r["start"]), int(r["end"])))

                        sequence = trim_positions(read.query_sequence, ranges)
                        click.echo(f"Processing read: {read.query_name}, ranges: {ranges}, original length: {read.query_length}, new_length: {len(sequence)}")
                        if len(sequence) > 0:
                            out.write(format_fasta_record(read.query_name, sequence))
                    else:
                        out.write(format_fasta_record(read.query_name, read.query_sequence))

    if not iter_exhausted(filters):
        sys.exit("WARNING: not all entries in the BED file were processed! Did you forget to sort them in the same order as the BAM file?")

cli.add_command(blastout_to_bed)
cli.add_command(filter_bam_to_fasta)

if __name__ == '__main__':
    cli()
