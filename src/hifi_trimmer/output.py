import bgzip
import csv
import polars as pl
import pysam
from hifi_trimmer.utils import trim_positions, format_fasta_record


def create_bed(actions: pl.LazyFrame, end_length: int) -> pl.LazyFrame:
    """Format the output from determine_actions() into valid BED."""
    return actions.select(
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


def write_region_file(hits: pl.LazyFrame) -> pl.LazyFrame:
    """Takes a pl.LazyFrame of determined adapter matches and returns
    a pl.LazyFrame with a single column containing the unique
    region specifications (adapter:start-end) for each hit.
    """
    return (
        hits.select(["sseqid", "sstart", "send"])
        .unique()
        .with_columns(
            pl.when(pl.col("sstart") > pl.col("send"))
            .then(
                pl.struct(
                    sstart="send",
                    send="sstart",
                )
            )
            .otherwise(pl.struct(["sstart", "send"]))
            .struct.field(["sstart", "send"]),
        )
        .select(
            pl.concat_str(
                pl.col("sseqid"),
                pl.lit(":"),
                pl.col("sstart"),
                pl.lit("-"),
                pl.col("send"),
            ).alias("region")
        )
    )


def filter_bam(outfile, bam, bed):
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
                    print(
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

    ## return the remaining filters - if we did all reads, should be empty
    return filters
