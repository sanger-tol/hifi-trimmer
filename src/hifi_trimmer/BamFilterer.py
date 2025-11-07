import bgzip
import click
import polars as pl
import pysam


class BamFilterer:
    def __init__(
        self, bam: click.File, bed: str, outfile: str, threads: int, fastq: bool
    ):
        self.threads = threads
        self.write_fastq = fastq

    def format_fastx_record(self, header: str, sequence: str, qual: str) -> str:
        """Format a header and sequence into FASTA format"""
        if qual is not None:
            return f"@{header}\n{sequence}\n+\n{qual}\n"
        else:
            return f">{header}\n{sequence}\n"

    def trim_positions(self, seq: str, ranges: list) -> str:
        """Trim DNA sequence seq to remove the positions specified in ranges.

        seq: string
        ranges: list of non-overlapping tuples with (start, end)
        """
        ## sort the ranges so we trim from the end - that way indexing stays constant
        ranges = sorted(ranges, key=lambda x: x[0], reverse=True)

        for start, end in ranges:
            seq = seq[:start] + seq[end:]

        return seq

    def filter_bam_with_bed(self, bam: click.File, bed: str, outfile: str):
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
            click.echo("WARN: BED file is empty! Reads will be streamed as-is.")
            filters = iter([])
            r = next(filters, None)

        with bgzip.BGZipWriter(outfile, num_threads=self.threads) as out:
            with pysam.AlignmentFile(
                bam, "rb", check_sq=False, require_index=False
            ) as b:
                ## Process: for each read in the BAM, check if it matches the current BED record.
                ## If yes, pull BED records until we reach a record for the next read.
                ## Then trim the sequence using the records pulled and write to fasta.
                ## If not, write record straight to disk
                for read in b.fetch(until_eof=True):
                    if r is not None and r["read"] == read.query_name:
                        ranges = [(int(r["start"]), int(r["end"]))]

                        while True:
                            r = next(filters, None)
                            if r is None or r["read"] != read.query_name:
                                break
                            ranges.append((int(r["start"]), int(r["end"])))

                        sequence = self.trim_positions(read.query_sequence, ranges)

                        qual = None
                        if self.write_fastq:
                            qual = self.trim_positions(read.qual, ranges)

                        print(
                            f"Processing read: {read.query_name}, ranges: {ranges}, original length: {read.query_length}, new_length: {len(sequence)}"
                        )
                        if len(sequence) > 0:
                            out.write(
                                self.format_fastx_record(
                                    read.query_name, sequence, qual
                                ).encode("utf-8")
                            )
                    else:
                        qual = None
                        if self.write_fastq:
                            qual = read.qual

                        out.write(
                            self.format_fastx_record(
                                read.query_name, read.query_sequence, qual
                            ).encode("utf-8")
                        )

        try:
            read = next(filters)
            print(f"Read {read['read']} not processed!")
            for read in filters:
                print(f"Read {read['read']} not processed!")

            raise RuntimeError(
                "ERROR: Not all entries in the BED file were processed! Are they sorted in the same order as the BAM file?"
            )
        except StopIteration:
            print("Read filtering complete!")
