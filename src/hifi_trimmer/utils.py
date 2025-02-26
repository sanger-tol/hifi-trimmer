import polars as pl


def check_records(blast: pl.LazyFrame, adapters: pl.DataFrame) -> bool:
    """Check if any adapters in a BLAST dataframe match to more than one
    adapter sequence in the YAML file.
    """
    counts = (
        blast.select(pl.col("sseqid"))
        .unique()
        .join_where(adapters, pl.col("sseqid").cast(pl.String).str.contains(pl.col("adapter")))
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
