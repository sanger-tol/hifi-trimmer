import polars as pl
import polars.selectors as cs


def match_hits(
    blastout: pl.LazyFrame, adapter_df: pl.LazyFrame, end_length: int
) -> pl.LazyFrame:
    """Determine whether each blast hit matches sufficient criteria for discard/trimming

    Keyword arguments:
    blastout: pl.DataFrame containing the results of a blast search with read lengths in column 13
    adapter_key: pl.DataFrame of adapters and their filtering/trimming parameters
    end_length: length of window at either end of read which is examined for trimming
    """
    ## When to discard a read due to adapter presence
    discard = (
        (pl.col("discard_middle")).and_(
            (pl.col("read_length") > 2 * end_length),
            (pl.col("qstart") > end_length),
            (pl.col("qend") < pl.col("read_length") - end_length),
            (pl.col("qend") > end_length),
            (pl.col("qstart") < pl.col("read_length") - end_length),
            (pl.col("pident") >= pl.col("middle_pident")),
            (pl.col("length") >= pl.col("middle_length")),
        )
    ) | (
        (pl.col("discard_end")).and_(
            (
                (pl.col("qend") < end_length)
                | (pl.col("qstart") > pl.col("read_length") - end_length)
            ),
            (pl.col("pident") >= pl.col("end_pident")),
            (pl.col("length") >= pl.col("end_length")),
        )
    )

    ## When to trim a read at the left end
    trim_l = (pl.col("trim_end")).and_(
        pl.col("qend") < end_length,
        (pl.col("pident") >= pl.col("end_pident")),
        (pl.col("length") >= pl.col("end_length")),
    )

    ## When to trim a read at the righr end
    trim_r = (pl.col("trim_end")).and_(
        pl.col("qstart") >= pl.col("read_length") - end_length,
        (pl.col("pident") >= pl.col("end_pident")),
        (pl.col("length") >= pl.col("end_length")),
    )

    return (
        blastout
        ## match each blast hit with an adapter
        .join_where(
            adapter_df, pl.col("sseqid").cast(pl.String).str.contains(pl.col("adapter"))
        )
        ## Determine which hits to keep and return rows with a "real" hit
        .with_columns(discard=discard, trim_l=trim_l, trim_r=trim_r)
        .filter(pl.any_horizontal("discard", "trim_l", "trim_r"))
    )


def determine_actions(
    hits: pl.LazyFrame, end_length: int, min_length: int
) -> pl.LazyFrame:
    """Take a dataframe of matches from match_hits()
    and process the result to decide which actions to take.

    Keyword arguments:
    hits: pl.DataFrame of blast hits with matches
    end_length: length of window at either end of read which is examined for trimming
    min_length: minimum length of a read after trimming to keep
    """
    hits_actions = hits.with_columns(
        (
            pl.col("read_length")
            - ((pl.col("trim_l").any() + pl.col("trim_r").any()) * end_length)
        )
        .over("qseqid")
        .alias("read_length_after_trimming"),
    ).with_columns(
        gen_action=(
            pl.when(
                pl.col("discard").any()
                | (pl.col("read_length_after_trimming") < min_length)
            )
            .then(pl.lit("discard"))
            .otherwise(pl.lit("trim"))
            .over("qseqid")
        )
    )

    discard = (
        hits_actions.filter(pl.col("gen_action") == "discard")
        .group_by("qseqid")
        .agg(
            pl.col("sseqid", "qstart", "qend", "read_length").get(
                pl.col("evalue").arg_min()
            ),
            action=pl.lit("discard"),
        )
    )

    trim_l = (
        hits_actions.filter(~(pl.col("gen_action") == "discard"), pl.col("trim_l"))
        .group_by("qseqid")
        .agg(
            pl.col("sseqid", "qstart", "qend", "read_length")
            .filter(pl.col("qend") <= end_length)
            .get(pl.col("evalue").filter(pl.col("qend") <= end_length).arg_min()),
            action=pl.lit("trim_l"),
        )
    )

    trim_r = (
        hits_actions.filter(~(pl.col("gen_action") == "discard"), pl.col("trim_r"))
        .group_by("qseqid")
        .agg(
            pl.col("sseqid", "qstart", "qend", "read_length")
            .filter(pl.col("qstart") >= pl.col("read_length") - end_length)
            .get(
                pl.col("evalue")
                .filter(pl.col("qstart") >= pl.col("read_length") - end_length)
                .arg_min()
            ),
            action=pl.lit("trim_r"),
        )
    )

    return pl.concat([discard, trim_l, trim_r])
