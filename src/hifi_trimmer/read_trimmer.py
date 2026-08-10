from pathlib import Path

import polars as pl
import pyfastx
import pysam
from loguru import logger

from hifi_trimmer.file_writers import BAMWriter, FastXWriter


class ReadTrimmer:
    def __init__(self, threads: int, in_format: str, out_format: str, out_bam_format_opts: list | None, gzip_fasta: bool = False):
        self.threads = threads
        self.in_format = in_format
        self.out_format = out_format

        # Setup SAM read mode
        if in_format in ["bam", "sam", "cram"]:
            if in_format == "bam":
                self.sam_read_mode = "rb"
            elif in_format == "cram":
                self.sam_read_mode = "rc"
            else:  # sam
                self.sam_read_mode = "r"
        else:
            self.sam_read_mode = None

        # Format options for output
        self.gzip_fasta = gzip_fasta
        self.bam_format_opts = (
            [opt.encode() for opt in out_bam_format_opts] if out_bam_format_opts else None
        )

    def _read_input(self, input_path):
        """Open input file based on format."""
        if self.in_format in ["sam", "bam", "cram"]:
            return pysam.AlignmentFile(
                input_path, self.sam_read_mode, check_sq=False, require_index=False
            )
        elif self.in_format == "fasta":
            return pyfastx.Fasta(str(input_path), build_index=False)
        elif self.in_format == "fastq":
            return pyfastx.Fastq(str(input_path), build_index=False)
        else:
            raise ValueError(f"Unknown input format: {self.in_format}")

    def _create_writer(self, output_file, input_header=None):
        """Create appropriate writer based on output format.

        Args:
            output_file: File-like object opened in binary mode (from click.File)
            input_header: SAM header dict (required for BAM/SAM/CRAM formats)
        """
        if self.out_format in ["sam", "bam", "cram"]:
            writer = BAMWriter(
                output_file,
                header=input_header,
                format_type=self.out_format,
                threads=self.threads,
                format_options=self.bam_format_opts,
            )
            writer.add_pg_line()
            return writer
        elif self.out_format in ["fasta", "fastq"]:
            return FastXWriter(output_file, format_type=self.out_format, gzip_output=self.gzip_fasta)
        else:
            raise ValueError(f"Unknown output format: {self.out_format}")

    def trim_positions(self, seq: str, ranges: list) -> str:
        """Trim a sequence (DNA or qual string) to remove the positions specified in ranges.

        Args:
            seq: String to trim
            ranges: List of non-overlapping tuples with (start, end)

        Returns:
            Trimmed sequence string
        """
        ranges = sorted(ranges, key=lambda x: x[0], reverse=True)
        for start, end in ranges:
            seq = seq[:start] + seq[end:]
        return seq

    def trim_input_with_bed(self, input_path: Path, bed: str, output_path):
        """
        Trim input sequences using a BED file specifying trim regions.

        Supports input/output in FASTA, FASTQ, SAM, BAM, or CRAM formats.

        Args:
            input_path: Path to input file
            bed: Path to BED file with columns: read, start, end, reason
            output_path: Output file path or file-like object
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

        input_file = self._read_input(input_path)
        input_header = None
        if self.in_format in ["sam", "bam", "cram"]:
            input_header = input_file.header.to_dict()

        writer = self._create_writer(output_path, input_header=input_header)

        if self.out_format in ["sam", "bam", "cram"]:
            with writer:
                self._process_pacbio_bam(input_file, writer, filters, r)
        else:  # FASTA/FASTQ
            if self.in_format in ["sam", "bam", "cram"]:
                self._process_bam_to_fastx(input_file, writer, filters, r)
            else:
                self._process_fastx_reads(input_file, writer, filters, r)

        # Verify all BED entries were processed
        try:
            read = next(filters)
            logger.warning(f"Read {read['read']} not processed!")
            for read in filters:
                logger.warning(f"Read {read['read']} not processed!")
            raise RuntimeError(
                "ERROR: Not all entries in the BED file were processed! "
                "Are they sorted in the same order as the input file?"
            )
        except StopIteration:
            logger.info("Read trimming complete!")

    def _process_pacbio_bam(self, input_file, writer, filters, r):
        """Process PacBio BAM (SAM/BAM/CRAM format)."""
        for read in input_file.fetch(until_eof=True):
            if r is not None and r["read"] == read.query_name:
                # Collect all ranges for this read
                ranges = [(int(r["start"]), int(r["end"]))]

                while True:
                    r = next(filters, None)
                    if r is None or r["read"] != read.query_name:
                        break
                    ranges.append((int(r["start"]), int(r["end"])))

                # Trim sequence and qualities
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
                    writer.write_record(read)
            else:
                writer.write_record(read)

    def _process_bam_to_fastx(self, input_file, writer, filters, r):
        """Convert BAM/SAM/CRAM to FASTA/FASTQ format."""
        for read in input_file.fetch(until_eof=True):
            name = read.query_name
            seq = read.query_sequence
            qual = read.qual if read.qual else None

            if r is not None and r["read"] == name:
                # Collect all ranges for this read
                ranges = [(int(r["start"]), int(r["end"]))]

                while True:
                    r = next(filters, None)
                    if r is None or r["read"] != name:
                        break
                    ranges.append((int(r["start"]), int(r["end"])))

                # Trim sequence and qualities
                new_seq = self.trim_positions(seq, ranges)
                new_qual = (
                    self.trim_positions(qual, ranges)
                    if qual
                    else None
                )

                logger.info(
                    f"Processing read: {name}, ranges: {ranges}, "
                    f"original length: {read.query_length}, new length: {len(new_seq)}",
                )

                if len(new_seq) > 0:
                    writer.write_record(name, new_seq, new_qual)
            else:
                # Write record as-is
                writer.write_record(name, seq, qual)

    def _process_fastx_reads(self, input_file, writer, filters, r):
        """Process FASTX reads (FASTA/FASTQ format)."""
        for record in input_file:
            # Unpack tuple from pyfastx (name, seq) for FASTA or (name, seq, qual) for FASTQ
            if self.in_format == "fasta":
                name, seq = record
                qual = None
            else:  # fastq
                name, seq, qual = record

            if r is not None and r["read"] == name:
                # Collect all ranges for this read
                ranges = [(int(r["start"]), int(r["end"]))]

                while True:
                    r = next(filters, None)
                    if r is None or r["read"] != name:
                        break
                    ranges.append((int(r["start"]), int(r["end"])))

                # Trim sequence and qualities (if FASTQ)
                new_seq = self.trim_positions(seq, ranges)
                new_qual = (
                    self.trim_positions(qual, ranges)
                    if qual
                    else None
                )

                logger.info(
                    f"Processing read: {name}, ranges: {ranges}, "
                    f"original length: {len(seq)}, new length: {len(new_seq)}",
                )

                if len(new_seq) > 0:
                    writer.write_record(name, new_seq, new_qual)
            else:
                # Write record as-is
                writer.write_record(name, seq, qual)
