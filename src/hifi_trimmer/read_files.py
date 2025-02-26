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
            "qseqid": pl.Categorical,
            "sseqid": pl.Categorical,
            "pident": pl.Float32,
            "length": pl.UInt8,
            "mismatch": pl.UInt8,
            "gapopen": pl.UInt8,
            "qstart": pl.UInt16,
            "qend": pl.UInt16,
            "sstart": pl.UInt8,
            "send": pl.UInt8,
            "evalue": pl.Float32,
            "bitscore": pl.Float32,
            "read_length": pl.UInt16,
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
