import polars as pl
import polars.selectors as cs


def summarise_blast(blast: pl.LazyFrame) -> pl.DataFrame:
    return (
        blast.group_by("sseqid")
        .len("n_hits")
        .rename({"sseqid": "adapter"})
        .sort("adapter")
        .collect()
    )


def summarise_hits(hits: pl.DataFrame) -> pl.DataFrame:
    return (
        hits.unpivot(
            cs.by_name(["discard", "trim_l", "trim_r"]),
            index=~cs.by_name(["discard", "trim_l", "trim_r"]),
            variable_name="action",
        )
        .filter(pl.col("value"))
        .group_by(["sseqid", "action"])
        .len(name="n_hits")
        .filter(pl.col("n_hits") > 0)
        .rename({"sseqid": "adapter"})
    )


def summarise_actions(actions: pl.DataFrame) -> pl.DataFrame:
    return (
        actions.with_columns(
            bases=pl.when(pl.col("action") == "discard")
            .then(pl.col("read_length"))
            .when(pl.col("action").str.contains("trim"))
            .then(pl.lit(150)),
            action=pl.when(pl.col("action") == "discard")
            .then(pl.lit("discard"))
            .otherwise(pl.lit("trim")),
        )
        .group_by(["sseqid", "action"])
        .agg(
            n_reads=pl.col("qseqid").unique().len(), bases_removed=pl.col("bases").sum()
        )
        .filter(pl.col("n_reads") > 0)
        .rename({"sseqid": "adapter"})
    )
