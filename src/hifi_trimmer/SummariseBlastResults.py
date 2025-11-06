import polars as pl
import polars.selectors as cs

from hifi_trimmer.BlastProcessor import BlastProcessor


class SummariseBlastResults:
    def __init__(self, result: BlastProcessor):
        self.blast = result.raw_blast
        self.hits = result.hits
        self.actions = result.actions
        self.end_length = result.end_length

    def _summarise_blast(self) -> pl.DataFrame:
        """
        Summarise a raw BLAST table, returning the number of hits by
        adapter.
        """
        return (
            self.blast.group_by("sseqid")
            .len("n_hits")
            .rename({"sseqid": "adapter"})
            .sort("adapter")
            .collect()
        )

    def _summarise_hits(self) -> pl.DataFrame:
        """
        Summarise a filtered BLAST table, returning the number of hits by
        adapter and action.
        """
        return (
            self.hits.select(
                sseqid=pl.col("sseqid"),
                discard=pl.col("discard"),
                trim=(pl.col("trim_l") + pl.col("trim_r")).cast(pl.Boolean),
            )
            .unpivot(
                cs.by_name(["discard", "trim"]),
                index="sseqid",
                variable_name="action",
            )
            .filter(pl.col("value"))
            .group_by(["sseqid", "action"])
            .len(name="n_hits")
            .filter(pl.col("n_hits") > 0)
            .rename({"sseqid": "adapter"})
            .sort(["adapter", "action"])
        )

    def _summarise_adapter_actions(self) -> pl.DataFrame:
        """
        Summarise an actions table, returning the number of reads affected
        and the number of bases removed for each adapter and action.
        """
        return (
            self.actions.with_columns(
                bases=pl.when(pl.col("action") == "discard")
                .then(pl.col("read_length"))
                .when(pl.col("action").str.contains("trim"))
                .then(pl.lit(self.end_length)),
                action=pl.when(pl.col("action") == "discard")
                .then(pl.lit("discard"))
                .otherwise(pl.lit("trim")),
            )
            .group_by(["sseqid", "action"])
            .agg(
                n_reads=pl.col("qseqid").unique().len(),
                bases_removed=pl.col("bases").sum(),
            )
            .filter(pl.col("n_reads") > 0)
            .rename({"sseqid": "adapter"})
            .sort(["adapter", "action"])
        )

    def summarise_actions(self) -> pl.DataFrame:
        """
        Summarise an actions table, returning the total number of reads affected
        and the number of bases removed for each action.
        """
        return (
            self.actions.with_columns(
                bases=pl.when(pl.col("action") == "discard")
                .then(pl.col("read_length"))
                .when(pl.col("action").str.contains("trim"))
                .then(pl.lit(self.end_length)),
                action=pl.when(pl.col("action") == "discard")
                .then(pl.lit("discard"))
                .otherwise(pl.lit("trim")),
            )
            .group_by(["action"])
            .agg(
                n_reads=pl.col("qseqid").unique().len(),
                bases_removed=pl.col("bases").sum(),
            )
        )

    def generate_summary(self) -> dict:
        """
        Takes the raw BLAST table, the BLAST table filtered for "valid" hits
        using criteria from the YAML, and the per-read "actions" table,
        and returns a dictionary with a summary about the adapter removal process:

        detected: how many hits per adapter detected in the raw BLAST output
        hits: for each adapter and action, the following stats:
            - number of blast hits
            - number of reads where this was the best hit
            - number of bases removed
        total_reads_discarded: total number of reads discarded
        total_reads_trimmed: total number of reads trimmed
        total_bases_removed: total length of removed bases
        """
        summary = {
            "detections": [],
            "hits": [],
            "total_reads_discarded": 0,
            "total_reads_trimmed": 0,
            "total_bases_removed": 0,
        }

        blast_summary = self._summarise_blast(self.blast)
        hits_summary = self._summarise_hits(self.hits)
        adapter_actions_summary = self._summarise_adapter_actions(self.actions)
        all_actions_summary = self._summarise_actions(self.actions)

        if not blast_summary.isEmpty():
            summary["detections"] = blast_summary.to_dicts()

            if not hits_summary.is_empty():
                summary["hits"] = (
                    hits_summary.join(adapter_actions_summary, on=["adapter", "action"])
                    .sort(["adapter", "action"])
                    .to_dicts()
                )

                summary["total_bases_removed"] = all_actions_summary[
                    "bases_removed"
                ].sum()
                summary["total_reads_discarded"] = all_actions_summary.filter(
                    pl.col("action") == "discard"
                )["n_reads"].sum()
                summary["total_reads_trimmed"] = all_actions_summary.filter(
                    pl.col("action") == "trim"
                )["n_reads"].sum()

        return summary
