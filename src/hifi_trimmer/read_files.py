import polars as pl
import yaml


def read_adapter_yaml(yaml_path: str) -> pl.LazyFrame:
    """Open an adapter YAML file and return it as a lazy pl.DataFrame"""
    try:
        adapters = yaml.safe_load(yaml_path)
    except yaml.YAMLError as exc:
        print(exc)

    return pl.DataFrame(
        adapters,
        schema={
            "adapter": pl.String,
            "discard_middle": pl.Boolean,
            "discard_end": pl.Boolean,
            "trim_end": pl.Boolean,
            "middle_pident": pl.Float32,
            "end_pident": pl.Float32,
            "middle_length": pl.UInt8,
            "end_length": pl.UInt8,
        },
    ).lazy()


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
            "pident": pl.Float32,
            "length": pl.UInt8,
            "mismatch": pl.UInt8,
            "gapopen": pl.UInt8,
            "qstart": pl.UInt32,
            "qend": pl.UInt32,
            "sstart": pl.UInt8,
            "send": pl.UInt8,
            "evalue": pl.Float32,
            "bitscore": pl.Float32,
            "read_length": pl.UInt32,
        },
    ).with_columns(
        pl.col("qseqid")
        .cast(pl.Categorical)
        .set_sorted(),  ## Set sorted for faster grouping operations
        pl.col("sseqid").cast(pl.Categorical),
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
