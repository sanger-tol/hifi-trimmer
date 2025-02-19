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


def read_readlengths(filtered_reads: set, bam: str) -> pl.DataFrame:
    """Count the lengths of reads in set filtered_reads using bamfile bam"""
    lengths = []
    with pysam.AlignmentFile(bam, "rb", check_sq=False, require_index=False) as b:
        for read in b.fetch(until_eof=True):
            if read.query_name in filtered_reads:
                lengths.append(
                    {"qseqid": read.query_name, "read_length": read.query_length}
                )

    return pl.DataFrame(lengths)


def read_blast(blast_path: str, bam: str = None) -> pl.LazyFrame:
    """Read BLAST outformat file to a LazyFrame.

    If a BAM file is provided and there is no 13th column (containing read lengths)
    then these can be quickly counted.

    Returns with the columns qseqid, sseqid, pident, length, qstart, qend, and evalue.
    qstart and qend are sorted to ignore information about strand.
    """
    blastout = pl.scan_csv(blast_path, has_header=False, separator="\t")

    ## if read lengths not provided but bam is - get them from the bam
    if len(blastout.head(1).collect().columns) == 12:
        if bam is not None:
            filtered_reads = set(
                blastout.select("column_1").unique().collect()["column_1"].to_list()
            )
            lengths = read_readlengths(filtered_reads, bam)
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
                "sstart",
                "send",
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
