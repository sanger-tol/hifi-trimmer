import click
import csv
import bgzip
from pathlib import Path
import polars as pl
import pysam
import yaml


def read_adapter_yaml(yaml_path: str) -> pl.LazyFrame:
    """Open an adapter YAML file and return it as a lazy pl.DataFrame"""
    try:
        adapters = yaml.safe_load(yaml_path)
    except yaml.YAMLError as exc:
        print(exc)

    return pl.DataFrame(adapters).lazy()


def read_blast(blast_path: str, bam: str = None) -> pl.LazyFrame:
    """Read BLAST outformat file to a LazyFrame.

    If a BAM file is provided and there is no 13th column (containing read lengths)
    then these can be quickly counted.

    Returns with the columns qseqid, sseqid, pident, length, qstart, qend, and evalue.
    qstart and qend are sorted to ignore information about strand.
    """
    blastout = pl.scan_csv(blast_path, has_header=False, separator="\t")

    ## if read lengths not provided but bam is - get them from the bam
    if len(blastout.head(1).collect().columns) == 12 and bam is not None:
        filtered_reads = set(
            blastout.select("column_1").unique().collect()["column_1"].to_list()
        )
        lengths = read_lengths(filtered_reads, bam)
        blastout = blastout.join(
            lengths.lazy(), left_on="column_1", right_on="qseqid", how="left", suffix=""
        )
    else:
        raise click.ClickException(
            "BLAST file only has 12 columns, but no BAM file has been provided to count!"
        )

    column_names = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "read_length",
    ]

    blastout = (
        blastout.rename(
            {
                old: new
                for old, new in zip(blastout.collect_schema().names(), column_names)
            }
        )
        ## reorder start-end to correct strandedness as we're not interested here
        .select(
            [
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "qstart",
                "qend",
                "evalue",
                "read_length",
            ]
        )
        .with_columns(
            pl.col("qseqid").set_sorted(),  ## Set sorted for faster grouping operations
            pl.when(pl.col("qstart") > pl.col("qend"))
            .then(
                pl.struct(
                    start="qend",
                    end="qstart",
                )
            )
            .otherwise(pl.struct(["qstart", "qend"]))
            .struct.field(["qstart", "qend"]),
        )
    )

    return blastout


def read_lengths(filtered_reads: set, bam: str) -> pl.DataFrame:
    """Count the lengths of reads in set filtered_reads using bamfile bam"""
    lengths = []
    with pysam.AlignmentFile(bam, "rb", check_sq=False, require_index=False) as b:
        for read in b.fetch(until_eof=True):
            if read.query_name in filtered_reads:
                lengths.append(
                    {"qseqid": read.query_name, "read_length": read.query_length}
                )

    return pl.DataFrame(lengths)


def check_records(blast: pl.LazyFrame, adapters: pl.DataFrame) -> bool:
    """Check if any adapters in a BLAST dataframe match to more than one
    adapter sequence in the YAML file.
    """
    counts = (
        blast.select(pl.col("sseqid"))
        .unique()
        .join_where(adapters, pl.col("sseqid").str.contains(pl.col("adapter")))
        .group_by("sseqid")
        .len()
        .filter(pl.col("len") > 1)
        .collect()["sseqid"]
        .to_list()
    )

    return counts


def format_fasta_record(header: str, sequence: str) -> str:
    """Format a header and sequence into FASTA format"""
    return f">{header}\n{sequence}\n"


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


def match_hits(
    blastout: pl.LazyFrame, adapter_df: pl.LazyFrame, end_length: int
) -> pl.LazyFrame:
    """Determine whether each blast hit matches sufficient criteria for discard/trimming

    Keyword arguments:
    blastout: pl.DataFrame containing the results of a blast search with read lengths in column 13
    adapter_key: pl.DataFrame of adapters and their filtering/trimming parameters
    end_length: length of window at either end of read which is examined for trimming
    """
    ## When to discard a read due to adapter presence
    discard = (
        (pl.col("discard_middle")).and_(
            (pl.col("read_length") > 2 * end_length),
            (pl.col("qstart") > end_length),
            (pl.col("qend") < pl.col("read_length") - end_length),
            (pl.col("pident") >= pl.col("middle_pident")),
            (pl.col("length") >= pl.col("middle_length")),
        )
    ) | (
        (pl.col("discard_end")).and_(
            (pl.col("qend") < end_length),
            (pl.col("qstart") > pl.col("read_length") - end_length),
            (pl.col("pident") >= pl.col("end_pident")),
            (pl.col("length") >= pl.col("end_length")),
        )
    )

    ## When to trim a read at the left end
    trim_l = (pl.col("trim_end")).and_(
        pl.col("qend") <= end_length,
        (pl.col("pident") >= pl.col("end_pident")),
        (pl.col("length") >= pl.col("end_length")),
    )

    ## When to trim a read at the righr end
    trim_r = (pl.col("trim_end")).and_(
        pl.col("qstart") >= pl.col("read_length") - end_length,
        (pl.col("pident") >= pl.col("end_pident")),
        (pl.col("length") >= pl.col("end_length")),
    )

    blastout = (
        blastout
        ## match each blast hit with an adapter
        ## below line more efficient but currently can't preserve order
        # .join_where(adapter_df, pl.col("sseqid").str.contains(pl.col("adapter")))
        .join(adapter_df, how="cross", maintain_order="left")
        .filter(pl.col("sseqid").str.contains(pl.col("adapter")))
        ## Determine which hits to keep and return rows with a "real" hit
        .with_columns(discard=discard, trim_l=trim_l, trim_r=trim_r)
        .filter(pl.any_horizontal("discard", "trim_l", "trim_r"))
    )

    return blastout


def determine_actions(
    blastout: pl.LazyFrame, end_length: int, min_length: int
) -> pl.LazyFrame:
    """Take a dataframe of matches from match_hits()
    and process the result to decide which actions to take.

    Keyword arguments:
    blastout: pl.DataFrame of blast hits with matches
    end_length: length of window at either end of read which is examined for trimming
    min_length: minimum length of a read after trimming to keep
    """

    def filter_get_min(getcol, maxcol, cond) -> pl.Expr:
        return getcol.filter(cond).get(maxcol.filter(cond).arg_min())

    ## Calculate read length after trimming (naiively) and then unpivot the trim columns so the
    ## trim section below operates separately for l and r
    result = (
        blastout.with_columns(
            (
                pl.col("read_length")
                - ((pl.col("trim_l").any() + pl.col("trim_r").any()) * end_length)
            )
            .over("qseqid")
            .alias("read_length_after_trimming"),
        )
        .with_columns(
            discard=(
                pl.col("discard") | (pl.col("read_length_after_trimming") < min_length)
            )
        )
        .drop("read_length_after_trimming")
        .group_by("qseqid", maintain_order=True)
        .agg(
            read_length=pl.col("read_length").first(),
            action=pl.when(pl.col("discard").any())
            .then(pl.lit(["discard"]))
            .when(pl.col("trim_r").any() & pl.col("trim_l").any())
            .then(pl.lit(["trim_l", "trim_r"]))
            .when(pl.col("trim_l").any() & ~pl.col("trim_r").any())
            .then(pl.lit(["trim_l"]))
            .when(pl.col("trim_r").any() & ~pl.col("trim_l").any())
            .then(pl.lit(["trim_r"]))
            .otherwise(pl.lit(["no_reason"])),
            cols=pl.when(pl.col("discard").any())
            .then(
                pl.concat_list(
                    pl.struct(["sseqid", "qstart", "qend"]).get(
                        pl.col("evalue").arg_max()
                    )
                )
            )
            .when(pl.col("trim_l").any() & pl.col("trim_r").any())
            .then(
                pl.concat_list(
                    filter_get_min(
                        pl.struct(["sseqid", "qstart", "qend"]),
                        pl.col("evalue"),
                        pl.col("qend") <= end_length,
                    ).implode(),
                    filter_get_min(
                        pl.struct(["sseqid", "qstart", "qend"]),
                        pl.col("evalue"),
                        pl.col("qstart") >= pl.col("read_length") - end_length,
                    ).implode(),
                )
            )
            .when(pl.col("trim_l").any() & ~pl.col("trim_r").any())
            .then(
                pl.concat_list(
                    filter_get_min(
                        pl.struct(["sseqid", "qstart", "qend"]),
                        pl.col("evalue"),
                        pl.col("qend") <= end_length,
                    )
                )
            )
            .when(pl.col("trim_r").any() & ~pl.col("trim_l").any())
            .then(
                pl.concat_list(
                    filter_get_min(
                        pl.struct(["sseqid", "qstart", "qend"]),
                        pl.col("evalue"),
                        pl.col("qstart") >= (pl.col("read_length") - end_length),
                    )
                )
            )
            .otherwise(
                pl.concat_list(
                    pl.struct(sseqid=pl.lit("none"), qstart=pl.lit(0), qend=pl.lit(0))
                )
            )
            .explode()  ## this could be tidied
            .explode(),
        )
        .explode("action", "cols")
        .with_columns(pl.col("cols").struct.field(["sseqid", "qstart", "qend"]))
    )

    return result


def create_bed(blastout: pl.LazyFrame, end_length: int) -> pl.LazyFrame:
    """Format the output from determine_actions() into valid BED."""
    return blastout.select(
        qseqid=pl.col("qseqid"),
        start=pl.when((pl.col("action") == "discard") | (pl.col("action") == "trim_l"))
        .then(pl.lit(0))
        .when(pl.col("action") == "trim_r")
        .then(pl.col("read_length") - end_length)
        .otherwise(pl.lit(0)),
        end=pl.when((pl.col("action") == "discard") | (pl.col("action") == "trim_r"))
        .then(pl.col("read_length"))
        .when(pl.col("action") == "trim_l")
        .then(pl.lit(end_length))
        .otherwise(pl.lit(0)),
        reason=pl.concat_str(pl.col("action"), pl.lit(":"), pl.col("sseqid")),
    )


@click.group()
def cli():
    """Main entry point for the tool."""
    pass


@click.command("blastout_to_bed")
@click.argument("blastout", type=click.File(mode="r"))
@click.argument("adapter_yaml", type=click.File(mode="r"))
@click.option(
    "-o",
    "--output",
    default="-",
    required=False,
    type=click.File(mode="wb"),
    help="Output file to write BED to.",
)
@click.option(
    "-b",
    "--bam",
    default=None,
    required=False,
    type=click.File(mode="rb"),
    help="If blastout file has no read length field, a BAM file of reads to get read lengths",
)
@click.option(
    "-ml",
    "--min_length_after_trimming",
    default=300,
    type=int,
    help="Minumum length of a read after trimming the ends in order not to be discarded.",
)
@click.option(
    "-el",
    "--end_length",
    default=150,
    type=int,
    help="Window size at either end of the read to be considered as 'ends' for searching.",
)
def blastout_to_bed(
    blastout, adapter_yaml, output, bam, min_length_after_trimming, end_length
):
    """Processes the input blastout file according to the adapter yaml key.

    BLASTOUT: tabular file resulting from BLAST with -outfmt "6 std qlen". If the qlen column
    is missing, lengths can be calculated by passing the --bam option.

    ADAPTER_YAML: yaml file contaning a list with the following fields per adapters:
 
        - name: (name of adapter. can be a regular expression) \n
        - discard_middle: True/False (discard read if adapter found in middle) \n
        - discard_end: True/False (discard read if adapter found in end) \n
        - trim_end: True/False (trim read if adapter found in end) \n
        - middle_pident: int (minimum pident requred to identify adapter in middle of read) \n
        - middle_length: int (minimum match length required to identify adapter in middle of read) \n
        - end_pident: int (minimum pident requred to identify adapter in end window) \n
        - end_length: int (minimum match length requred to identify adapter in end window) \n

    Output:
    By default, writes BED to standard output. This can be redirected with the -o/--output option.
    """
    adapters = read_adapter_yaml(adapter_yaml)
    blast = read_blast(blastout, bam)

    ## Make sure each barcode that matches an adapter matches only one adapter
    check = check_records(blast, adapters)
    if len(check) > 0:
        raise click.ClickException(
            f"Barcodes {', '.join(check)} match to more than one adapter in {adapter_yaml.name}!"
        )

    ## process the blastout file
    hits = match_hits(blast, adapters, end_length)
    actions = determine_actions(hits, end_length, min_length_after_trimming)
    bed = create_bed(actions, end_length)

    bed.collect().write_csv(output, separator="\t", include_header=False)


@click.command("filter_bam_to_fasta")
@click.argument("bed", type=click.File(mode="r"))
@click.argument("bam", type=click.File(mode="rb"))
@click.argument("outfile", default=None, required=True, type=click.File(mode="wb"))
def filter_bam_to_fasta(bed, bam, outfile):
    """
    Filter the reads stored in a BAM file using the appropriate BED file produced
    by blastout_to_bed and write to a bgzipped fasta file.

    BED: BED file describing regions of the read set to exclude.
    BAM: BAM file in which to filter reads
    OUTFILE: File to write the filtered reads to (bgzipped).
    """
    if outfile is None:
        outfile = Path(bam).stem + ".filtered.fa.gz"

    ## read BED as dictionary
    filters = csv.DictReader(
        bed, delimiter="\t", fieldnames=["read", "start", "end", "reason"]
    )

    with bgzip.BGZipWriter(outfile) as out:
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
                    click.echo(
                        f"Processing read: {read.query_name}, ranges: {ranges}, original length: {read.query_length}, new_length: {len(sequence)}"
                    )
                    if len(sequence) > 0:
                        out.write(
                            format_fasta_record(read.query_name, sequence).encode(
                                "utf-8"
                            )
                        )
                else:
                    out.write(
                        format_fasta_record(
                            read.query_name, read.query_sequence
                        ).encode("utf-8")
                    )

    try:
        read = next(filters)
        print(f"Read {read['read']} not processed!")
        for read in filters:
            print(f"Read {read['read']} not processed!")

        raise click.ClickException(
            "Not all entries in the BED file were processed! Did you forget to sort them in the same order as the BAM file?"
        )
    except StopIteration:
        print("Read filtering complete!")


cli.add_command(blastout_to_bed)
cli.add_command(filter_bam_to_fasta)

if __name__ == "__main__":
    cli()
