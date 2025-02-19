import polars as pl


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
            (pl.col("qend") > end_length),
            (pl.col("qstart") < pl.col("read_length") - end_length),
            (pl.col("pident") >= pl.col("middle_pident")),
            (pl.col("length") >= pl.col("middle_length")),
        )
    ) | (
        (pl.col("discard_end")).and_(
            (pl.col("qend") < end_length),
            (pl.col("qstart") > pl.col("read_length") - end_length),
            (pl.col("pident") >= pl.col("end_pident")),
            (pl.col("length") >= pl.col("end_length")),
        )
    )

    ## When to trim a read at the left end
    trim_l = (pl.col("trim_end")).and_(
        pl.col("qend") <= end_length,
        (pl.col("pident") >= pl.col("end_pident")),
        (pl.col("length") >= pl.col("end_length")),
    )

    ## When to trim a read at the righr end
    trim_r = (pl.col("trim_end")).and_(
        pl.col("qstart") >= pl.col("read_length") - end_length,
        (pl.col("pident") >= pl.col("end_pident")),
        (pl.col("length") >= pl.col("end_length")),
    )

    blastout = (
        blastout
        ## match each blast hit with an adapter
        ## below line more efficient but currently can't preserve order
        # .join_where(adapter_df, pl.col("sseqid").str.contains(pl.col("adapter")))
        .join(adapter_df, how="cross", maintain_order="left")
        .filter(pl.col("sseqid").str.contains(pl.col("adapter")))
        ## Determine which hits to keep and return rows with a "real" hit
        .with_columns(discard=discard, trim_l=trim_l, trim_r=trim_r)
        .filter(pl.any_horizontal("discard", "trim_l", "trim_r"))
    )

    return blastout


def determine_actions(
    blastout: pl.LazyFrame, end_length: int, min_length: int
) -> pl.LazyFrame:
    """Take a dataframe of matches from match_hits()
    and process the result to decide which actions to take.

    Keyword arguments:
    blastout: pl.DataFrame of blast hits with matches
    end_length: length of window at either end of read which is examined for trimming
    min_length: minimum length of a read after trimming to keep
    """

    def filter_get_min(getcol, maxcol, cond) -> pl.Expr:
        return getcol.filter(cond).get(maxcol.filter(cond).arg_min())

    ## Calculate read length after trimming (naiively) and then unpivot the trim columns so the
    ## trim section below operates separately for l and r
    result = (
        blastout.with_columns(
            (
                pl.col("read_length")
                - ((pl.col("trim_l").any() + pl.col("trim_r").any()) * end_length)
            )
            .over("qseqid")
            .alias("read_length_after_trimming"),
        )
        .with_columns(
            discard=(
                pl.col("discard") | (pl.col("read_length_after_trimming") < min_length)
            )
        )
        .drop("read_length_after_trimming")
        .group_by("qseqid", maintain_order=True)
        .agg(
            read_length=pl.col("read_length").first(),
            action=pl.when(pl.col("discard").any())
            .then(pl.lit(["discard"]))
            .when(pl.col("trim_r").any() & pl.col("trim_l").any())
            .then(pl.lit(["trim_l", "trim_r"]))
            .when(pl.col("trim_l").any() & ~pl.col("trim_r").any())
            .then(pl.lit(["trim_l"]))
            .when(pl.col("trim_r").any() & ~pl.col("trim_l").any())
            .then(pl.lit(["trim_r"]))
            .otherwise(pl.lit(["no_reason"])),
            cols=pl.when(pl.col("discard").any())
            .then(
                pl.concat_list(
                    pl.struct(["sseqid", "qstart", "qend"]).get(
                        pl.col("evalue").arg_max()
                    )
                )
            )
            .when(pl.col("trim_l").any() & pl.col("trim_r").any())
            .then(
                pl.concat_list(
                    filter_get_min(
                        pl.struct(["sseqid", "qstart", "qend"]),
                        pl.col("evalue"),
                        pl.col("qend") <= end_length,
                    ).implode(),
                    filter_get_min(
                        pl.struct(["sseqid", "qstart", "qend"]),
                        pl.col("evalue"),
                        pl.col("qstart") >= pl.col("read_length") - end_length,
                    ).implode(),
                )
            )
            .when(pl.col("trim_l").any() & ~pl.col("trim_r").any())
            .then(
                pl.concat_list(
                    filter_get_min(
                        pl.struct(["sseqid", "qstart", "qend"]),
                        pl.col("evalue"),
                        pl.col("qend") <= end_length,
                    )
                )
            )
            .when(pl.col("trim_r").any() & ~pl.col("trim_l").any())
            .then(
                pl.concat_list(
                    filter_get_min(
                        pl.struct(["sseqid", "qstart", "qend"]),
                        pl.col("evalue"),
                        pl.col("qstart") >= (pl.col("read_length") - end_length),
                    )
                )
            )
            .otherwise(
                pl.concat_list(
                    pl.struct(sseqid=pl.lit("none"), qstart=pl.lit(0), qend=pl.lit(0))
                )
            )
            .explode()  ## this could be tidied
            .explode(),
        )
        .explode("action", "cols")
        .with_columns(pl.col("cols").struct.field(["sseqid", "qstart", "qend"]))
    )

    return result
