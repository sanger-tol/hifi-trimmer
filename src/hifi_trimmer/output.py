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
    blast: pl.LazyFrame, hits: pl.LazyFrame, bed: pl.LazyFrame
) -> pl.LazyFrame:
    blast_summary = (
        blast.group_by("sseqid")
        .len("n_hits")
        .rename({"sseqid": "adapter"})
        .collect()
        .to_dict(as_series=False)
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
        .sort(["adapter", "action"])
        .collect()
        .to_dict(as_series=False)
    )

    n_bases_removed = (
        bed.select(n_bases=pl.col("end") - pl.col("start")).collect()["n_bases"].sum()
    )

    return {
        "detections": blast_summary,
        "hits": hit_summary,
        "bases_removed": n_bases_removed,
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
