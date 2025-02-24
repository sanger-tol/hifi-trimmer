import click
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


def read_blast(blast_path: str) -> pl.LazyFrame:
    """Read BLAST outformat file to a LazyFrame.

    If a BAM file is provided and there is no 13th column (containing read lengths)
    then these can be quickly counted.

    Returns with the columns qseqid, sseqid, pident, length, qstart, qend, and evalue.
    qstart and qend are sorted to ignore information about strand.
    """
    blastout = pl.scan_csv(
        blast_path,
        has_header=False,
        separator="\t",
        schema={
            "qseqid": pl.String,
            "sseqid": pl.String,
            "pident": pl.Float64,
            "length": pl.Int64,
            "mismatch": pl.Int64,
            "gapopen": pl.Int64,
            "qstart": pl.Int64,
            "qend": pl.Int64,
            "sstart": pl.Int64,
            "send": pl.Int64,
            "evalue": pl.Float64,
            "bitscore": pl.Float64,
            "read_length": pl.Int64,
        },
    ).with_columns(
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

    return blastout
