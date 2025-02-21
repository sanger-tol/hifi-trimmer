import bgzip
import csv
import polars as pl
import polars.selectors as cs
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


def write_summary(
    blast: pl.LazyFrame, hits: pl.LazyFrame, actions: pl.LazyFrame
) -> pl.LazyFrame:
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
    blast_summary = (
        blast.group_by("sseqid")
        .len("n_hits")
        .rename({"sseqid": "adapter"})
        .sort("adapter")
        .collect()
        .to_dicts()
    )

    hit_summary = (
        hits.unpivot(
            cs.by_name(["discard", "trim_l", "trim_r"]),
            index=~cs.by_name(["discard", "trim_l", "trim_r"]),
            variable_name="action",
        )
        .filter(pl.col("value"))
        .group_by(["sseqid", "action"])
        .len(name="n_hits")
        .rename({"sseqid": "adapter"})
    )

    actions_summary = (
        actions.with_columns(
            bases=pl.when(pl.col("action").str.contains("discard"))
            .then(pl.col("read_length").first())
            .when(pl.col("action").str.contains("trim"))
            .then(pl.lit(150))
            .otherwise(pl.lit(0))
        )
        .group_by(["sseqid", "action"])
        .agg(n_reads = pl.col("sseqid").len(),
             bases_removed=pl.col("bases").sum())
        .rename({"sseqid": "adapter"})
    )

    hit_actions_summary = (
        hit_summary.join(actions_summary, how="left", on=["adapter", "action"])
        .with_columns(
            n_reads=pl.col("n_reads").fill_null(strategy="zero"),
            bases_removed=pl.col("bases_removed").fill_null(strategy="zero"))
        .sort(["adapter", "action"])
        .collect()
    )

    n_bases_removed = hit_actions_summary["bases_removed"].sum()
    n_reads_discarded = (hit_actions_summary
        .filter(pl.col("action") == "discard")
        ["n_reads"].sum()
    )

    n_reads_trimmed = (actions
        .filter(pl.col("action").str.contains("trim"))
        .select(len = pl.col("qseqid").unique().len())
        .collect()["len"].to_list()[0]
    )

    print(n_reads_trimmed)

    return {
        "detections": blast_summary,
        "hits": hit_actions_summary.to_dicts(),
        "total_reads_discarded": n_reads_discarded,
        "total_reads_trimmed": n_reads_trimmed,
        "total_bases_removed": n_bases_removed,
    }


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
    )


def filter_bam_with_bed(outfile, bam, bed, threads):
    ## read BED as dictionary

    with bgzip.BGZipReader(bed, num_threads=threads) as file:
        filters = csv.DictReader(
            file.read().tobytes().decode("utf-8").splitlines(),
            delimiter="\t",
            fieldnames=["read", "start", "end", "reason"],
        )

    with bgzip.BGZipWriter(outfile, num_threads=threads) as out:
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
