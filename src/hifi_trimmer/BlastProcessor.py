from pathlib import Path

import click
import polars as pl
import yaml


class BlastProcessor:
    def __init__(self, blast, yaml, end_length=150, min_length=300):
        ## Initialise params
        self.blast_file = Path(blast).resolve()
        self.adapter_yaml_file = Path(yaml).resolve()
        self.end_length = end_length
        self.min_length = min_length

        ## Initialise default empty data structures
        self.adapter_yaml = None
        self.raw_blast_table = None
        self.hits = None
        self.actions = None
        self.bed = b""

        ## Check if BLAST file is empty
        try:
            self._read_blast(blast).head(1).collect()
        except pl.exceptions.NoDataError:
            click.echo(f"WARNING: {blast} was empty! BED output will be empty.")
            return

        ## Read raw data
        self.raw_blast_table = self._read_blast(blast)
        self.adapter_table = self._read_adapter_yaml(yaml)

        ## Validate inputs
        self._check_duplicate_adapter_matches(self.raw_blast_table, self.adapter_table)

        ## Determine hits
        self.hits = self._determine_blast_hits(
            self.raw_blast_table, self.adapter_table
        ).collect()

        if self.hits.is_empty():
            self.hits = None
            click.echo(
                "WARNING: After filtering, no BLAST hits remain. BED output will be empty."
            )
            return

        # Generate actions and BED output
        self.actions = self._determine_per_read_actions(self.hits.lazy()).collect()
        self.bed = (
            self.generate_bed(self.actions.lazy())
            .collect()
            .write_csv(separator="\t", include_header=False)
            .encode("utf-8")
        )

    def _read_adapter_yaml(self, yaml_path: str) -> pl.LazyFrame:
        """Open an adapter YAML file and return it as a lazy pl.DataFrame"""
        try:
            with open(yaml_path) as f:
                adapters = yaml.safe_load(f)
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

    def _read_blast(self, blast_path: str) -> pl.LazyFrame:
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

    def _check_duplicate_adapter_matches(
        self, blast: pl.LazyFrame, adapter_df: pl.DataFrame
    ) -> list:
        """Check if any adapters in a BLAST dataframe match to more than one
        adapter sequence in the YAML file.
        """
        counts = (
            blast.select(pl.col("sseqid"))
            .unique()
            .with_columns(pl.col("sseqid").cast(pl.String))
            .join_where(adapter_df, pl.col("sseqid").str.contains(pl.col("adapter")))
            .group_by("sseqid")
            .len()
            .filter(pl.col("len") > 1)
            .collect()["sseqid"]
            .to_list()
        )

        if len(counts) > 0:
            raise ValueError(
                f"Adapters {', '.join(counts)} in the BLAST table match to more than one adapter defined in the YAML!"
            )

    def _determine_blast_hits(
        self, blast: pl.LazyFrame, adapter_table: pl.LazyFrame
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
                (pl.col("read_length") > 2 * self.end_length),
                (pl.col("qstart") > self.end_length),
                (pl.col("qend") < pl.col("read_length") - self.end_length),
                (pl.col("qend") > self.end_length),
                (pl.col("qstart") < pl.col("read_length") - self.end_length),
                (pl.col("pident") >= pl.col("middle_pident")),
                (pl.col("length") >= pl.col("middle_length")),
            )
        ) | (
            (pl.col("discard_end")).and_(
                (
                    (pl.col("qend") < self.end_length)
                    | (pl.col("qstart") > pl.col("read_length") - self.end_length)
                ),
                (pl.col("pident") >= pl.col("end_pident")),
                (pl.col("length") >= pl.col("end_length")),
            )
        )

        ## When to trim a read at the left end
        trim_l = (pl.col("trim_end")).and_(
            pl.col("qend") < self.end_length,
            (pl.col("pident") >= pl.col("end_pident")),
            (pl.col("length") >= pl.col("end_length")),
        )

        ## When to trim a read at the righr end
        trim_r = (pl.col("trim_end")).and_(
            pl.col("qstart") >= pl.col("read_length") - self.end_length,
            (pl.col("pident") >= pl.col("end_pident")),
            (pl.col("length") >= pl.col("end_length")),
        )

        return (
            blast
            ## match each blast hit with an adapter
            .join_where(
                adapter_table,
                pl.col("sseqid").cast(pl.String).str.contains(pl.col("adapter")),
            )
            ## Determine which hits to keep and return rows with a "real" hit
            .with_columns(discard=discard, trim_l=trim_l, trim_r=trim_r)
            .filter(pl.any_horizontal("discard", "trim_l", "trim_r"))
        )

    def _determine_per_read_actions(self, hits: pl.LazyFrame) -> pl.LazyFrame:
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
                - ((pl.col("trim_l").any() + pl.col("trim_r").any()) * self.end_length)
            )
            .over("qseqid")
            .alias("read_length_after_trimming"),
        ).with_columns(
            gen_action=(
                pl.when(
                    pl.col("discard").any()
                    | (pl.col("read_length_after_trimming") < self.min_length)
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
                .filter(pl.col("qend") <= self.end_length)
                .get(
                    pl.col("evalue").filter(pl.col("qend") <= self.end_length).arg_min()
                ),
                action=pl.lit("trim_l"),
            )
        )

        trim_r = (
            hits_actions.filter(~(pl.col("gen_action") == "discard"), pl.col("trim_r"))
            .group_by("qseqid")
            .agg(
                pl.col("sseqid", "qstart", "qend", "read_length")
                .filter(pl.col("qstart") >= pl.col("read_length") - self.end_length)
                .get(
                    pl.col("evalue")
                    .filter(pl.col("qstart") >= pl.col("read_length") - self.end_length)
                    .arg_min()
                ),
                action=pl.lit("trim_r"),
            )
        )

        return pl.concat([discard, trim_l, trim_r])

    def generate_bed(self, actions) -> pl.LazyFrame:
        """Format the output from determine_actions() into valid BED."""
        return actions.select(
            qseqid=pl.col("qseqid"),
            start=pl.when(
                (pl.col("action") == "discard") | (pl.col("action") == "trim_l")
            )
            .then(pl.lit(0))
            .when(pl.col("action") == "trim_r")
            .then(pl.col("read_length") - self.end_length)
            .otherwise(pl.lit(0)),
            end=pl.when(
                (pl.col("action") == "discard") | (pl.col("action") == "trim_r")
            )
            .then(pl.col("read_length"))
            .when(pl.col("action") == "trim_l")
            .then(pl.lit(self.end_length))
            .otherwise(pl.lit(0)),
            reason=pl.concat_str(pl.col("action"), pl.lit(":"), pl.col("sseqid")),
        ).sort("qseqid", "start")

    def generate_hits(self) -> pl.LazyFrame:
        return self.hits.select(
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
