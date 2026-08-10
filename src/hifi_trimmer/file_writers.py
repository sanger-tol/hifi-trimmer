import bgzip
import sys
from importlib.metadata import version
from typing import IO

import pysam


class FastXWriter:
    def __init__(self, output_file: IO[bytes], format_type: str, gzip_output: bool = False):
        """
        Args:
            output_file: File-like object opened in binary mode (from click.File)
            format_type: "fasta" or "fastq"
        """
        self.output_file = output_file
        self.format_type = format_type

        if gzip_output:
            self.output_file = bgzip.BGZipWriter(fileobj=output_file)


    def write_record(self, name: str, sequence: str, quality: str | None = None):
        """Write a FASTA or FASTQ record."""
        if self.format_type == "fasta":
            self._write_fasta_record(name, sequence)
        else:  # fastq
            self._write_fastq_record(name, sequence, quality)

    def _write_fasta_record(self, name: str, sequence: str):
        """Write a FASTA record."""
        self.output_file.write(f">{name}\n{sequence}\n".encode())

    def _write_fastq_record(self, name: str, sequence: str, quality: str | None = None):
        """Write a FASTQ record."""
        if quality is None:
            quality = "!" * len(sequence)
        self.output_file.write(f"@{name}\n{sequence}\n+\n{quality}\n".encode())


class BAMWriter:
    """Writer for BAM/SAM/CRAM formats with trimming support and PG header management."""

    def __init__(
        self,
        output_path,
        header: dict | None,
        format_type: str,
        threads: int = 1,
        format_options: list | None= None,
    ):
        """
        Args:
            output_path: Output file path or file-like object
            header: SAM header dict
            format_type: "sam", "bam", or "cram"
            threads: Number of threads for BAM compression
            format_options: Format-specific options (encoded as bytes)
        """
        self.output_path = output_path

        if header is None:
            raise ValueError("Input file has no header!")

        self.header = header
        self.format_type = format_type
        self.threads = threads
        self.format_options = format_options
        self.file_handle = None

    def _get_output_mode(self) -> str:
        """Get pysam output mode based on format type."""
        if self.format_type == "bam":
            return "wb"
        elif self.format_type == "cram":
            return "wc"
        else:  # sam
            return "w"

    def add_pg_line(self):
        """Add a PG line for hifi-trimmer to the header."""
        pg_entry = {
            "ID": "hifi-trimmer",
            "PN": "hifi-trimmer",
            "CL": " ".join(sys.argv),
            "VN": version("hifi-trimmer"),
        }

        if "PG" in self.header and len(self.header["PG"]) > 0:
            pg_entry["PP"] = self.header["PG"][-1]["ID"]

        self.header.setdefault("PG", []).append(pg_entry)

    def open(self):
        """Open the output file for writing."""
        output_mode = self._get_output_mode()
        self.file_handle = pysam.AlignmentFile(
            self.output_path,
            output_mode,
            header=pysam.AlignmentHeader.from_dict(self.header),
            threads=self.threads,
            format_options=self.format_options,
        )

    def close(self):
        """Close the output file."""
        if self.file_handle:
            self.file_handle.close()

    def write_record(self, aligned_segment):
        """Write a BAM record."""
        if self.file_handle:
            self.file_handle.write(aligned_segment)

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
