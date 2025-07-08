import bgzip
import click
import polars as pl
import pysam

from hifi_trimmer.utils import trim_positions, format_fastx_record


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
    ).sort("qseqid", "start")


def write_summary(
    blast_summary: pl.DataFrame,
    hits_summary: pl.DataFrame,
    actions_summary: pl.LazyFrame,
) -> dict:
    """
    Takes the raw BLAST table, the BLAST table filtered for "valid" hits
    using criteria from the YAML, and the per-read "actions" table,
    and returns a dictionary with a summary about the adapter removal process:

    detected: how many hits per adapter detected in the raw BLAST output
    hits: for each adapter and action, the following stats:
        - number of blast hits
        - number of reads where this was the best hit
        - number of bases removed
    total_reads_discarded: total number of reads discarded
    total_reads_trimmed: total number of reads trimmed
    total_bases_removed: total length of removed bases
    """
    summary = {
        "detections": [],
        "hits": [],
        "total_reads_discarded": 0,
        "total_reads_trimmed": 0,
        "total_bases_removed": 0,
    }

    if blast_summary is not None:
        summary["detections"] = blast_summary.to_dicts()

        if hits_summary is not None and not hits_summary.is_empty():
            summary["hits"] = (
                hits_summary.join(actions_summary, on=["adapter", "action"])
                .sort(["adapter", "action"])
                .to_dicts()
            )

            summary["total_bases_removed"] = actions_summary["bases_removed"].sum()
            summary["total_reads_discarded"] = actions_summary.filter(
                pl.col("action") == "discard"
            )["n_reads"].sum()
            summary["total_reads_trimmed"] = actions_summary.filter(
                pl.col("action") == "trim"
            )["n_reads"].sum()

    return summary


def write_hits(hits: pl.LazyFrame) -> pl.LazyFrame:
    return hits.select(
        [
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
            "discard",
            "trim_l",
            "trim_r",
        ]
    ).collect()


def filter_bam_with_bed(
    outfile: str, bam: click.File, bed: str, threads: int, fastq: bool
) -> iter:
    try:
        bed_df = pl.read_csv(
            bed,
            separator="\t",
            has_header=False,
            schema={
                "read": pl.String,
                "start": pl.Int64,
                "end": pl.Int64,
                "reason": pl.String,
            },
        )
        filters = bed_df.iter_rows(named=True)

        r = next(filters, None)
    except pl.exceptions.NoDataError:
        click.echo("WARN: BED file is empty! Reads will be streamed as-is.")
        filters = iter([])
        r = next(filters, None)

    with bgzip.BGZipWriter(outfile, num_threads=threads) as out:
        with pysam.AlignmentFile(bam, "rb", check_sq=False, require_index=False) as b:
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

                    qual = None
                    if fastq:
                        qual = trim_positions(read.qual, ranges)

                    print(
                        f"Processing read: {read.query_name}, ranges: {ranges}, original length: {read.query_length}, new_length: {len(sequence)}"
                    )
                    if len(sequence) > 0:
                        out.write(
                            format_fastx_record(read.query_name, sequence, qual).encode(
                                "utf-8"
                            )
                        )
                else:
                    qual = None
                    if fastq:
                        qual = read.qual

                    out.write(
                        format_fastx_record(
                            read.query_name, read.query_sequence, qual
                        ).encode("utf-8")
                    )

    ## return the remaining filters - if we did all reads, should be empty
    return filters
