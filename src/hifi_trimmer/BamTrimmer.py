import sys
from importlib.metadata import version
from pathlib import Path

import click
import polars as pl
import pysam
from loguru import logger


class BamTrimmer:
    def __init__(self, threads: int, format: str, format_opts: list | None):
        self.threads = threads
        self.output_mode = "w"
        if format == "bam":
            self.output_mode = "wb"
        elif format == "cram":
            self.output_mode = "wc"
        self.format_opts = (
            [opt.encode() for opt in format_opts] if format_opts else None
        )

    def trim_positions(self, seq: str, ranges: list) -> str:
        """Trim a sequence (DNA or qual string) to remove the positions specified in ranges.

        seq: string
        ranges: list of non-overlapping tuples with (start, end)
        """
        ranges = sorted(ranges, key=lambda x: x[0], reverse=True)
        for start, end in ranges:
            seq = seq[:start] + seq[end:]
        return seq

    def generate_PG_line(self, header: dict):
        """Add a PG line for hifi-trimmer to a BAM header dict."""

        pg_entry = {
            "ID": "hifi-trimmer",
            "PN": "hifi-trimmer",
            "CL": " ".join(sys.argv),
            "VN": version("hifi-trimmer"),
        }

        if "PG" in header and len(header["PG"]) > 0:
            pg_entry["PP"] = header["PG"][-1]["ID"]

        header.setdefault("PG", []).append(pg_entry)

        return header

    def trim_bam_with_bed(self, bam: click.File, bed: str, outfile):
        """
        Filter/trim an unaligned BAM file using a BED file and write BAM output.

        Args:
            bam: Input BAM file handle (from click).
            bed: Path to BED file with columns: read, start, end, reason.
            outfile: Output file handle (e.g. sys.stdout.buffer or a click.File('wb')).
                     Passed directly to pysam.AlignmentFile.
        """
        try:
            bed_df = pl.read_csv(
                bed,
                separator="\t",
                has_header=False,
                schema={
                    "read": pl.String,
                    "start": pl.Int64,
                    "end": pl.Int64,
                    "reason": pl.String,
                },
            )
            filters = bed_df.iter_rows(named=True)
            r = next(filters, None)
        except pl.exceptions.NoDataError:
            logger.warning("WARN: BED file is empty! Reads will be streamed as-is.")
            filters = iter([])
            r = next(filters, None)

        extension = Path(bam.name).suffix[1:].lower()
        in_read_mode = "r"
        if extension == "bam":
            in_read_mode = "rb"
        elif extension == "cram":
            in_read_mode = "rc"

        with pysam.AlignmentFile(
            bam, in_read_mode, check_sq=False, require_index=False
        ) as b:
            header = self.generate_PG_line(b.header.to_dict())

            with pysam.AlignmentFile(
                outfile,
                self.output_mode,
                header=pysam.AlignmentHeader.from_dict(header),
                threads=self.threads,
                format_options=self.format_opts,
            ) as out:
                for read in b.fetch(until_eof=True):
                    if r is not None and r["read"] == read.query_name:
                        ranges = [(int(r["start"]), int(r["end"]))]

                        while True:
                            r = next(filters, None)
                            if r is None or r["read"] != read.query_name:
                                break
                            ranges.append((int(r["start"]), int(r["end"])))

                        new_seq = self.trim_positions(read.query_sequence, ranges)
                        new_qual = (
                            self.trim_positions(read.qual, ranges)
                            if read.qual
                            else None
                        )

                        logger.info(
                            f"Processing read: {read.query_name}, ranges: {ranges}, "
                            f"original length: {read.query_length}, new length: {len(new_seq)}",
                        )

                        if len(new_seq) > 0:
                            read.query_sequence = new_seq
                            read.query_qualities = (
                                pysam.qualitystring_to_array(new_qual)
                                if new_qual
                                else None
                            )
                            out.write(read)
                    else:
                        out.write(read)

        try:
            read = next(filters)
            logger.warning(f"Read {read['read']} not processed!")
            for read in filters:
                logger.warning(f"Read {read['read']} not processed!")
            raise RuntimeError(
                "ERROR: Not all entries in the BED file were processed! "
                "Are they sorted in the same order as the BAM file?"
            )
        except StopIteration:
            logger.info("Read filtering complete!")
