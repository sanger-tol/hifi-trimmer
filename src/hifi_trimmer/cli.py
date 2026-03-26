import sys
from pathlib import Path

import bgzip
import click
from loguru import logger

from hifi_trimmer.BamTrimmer import BamTrimmer
from hifi_trimmer.BlastProcessor import BlastProcessor
from hifi_trimmer.SummariseBlastResults import SummariseBlastResults


@click.group()
@click.version_option(package_name="hifi_trimmer")
def cli():
    """Main entry point for the tool."""
    logger.remove()
    logger.add(
        sys.stderr,
        format="[{time:HH:mm:ss}] | hifi-trimmer | {level} - {message}",
        level="INFO",
    )


@click.command("process-blast")
@click.argument("blastout", type=click.Path(exists=True))
@click.argument("adapter-yaml", type=click.Path(exists=True))
@click.option(
    "-p",
    "--prefix",
    default=None,
    required=False,
    type=str,
    help="Output prefix for results. Defaults to the basename of the blastout if not provided.",
)
@click.option(
    "-ml",
    "--min-length-after-trimming",
    default=300,
    show_default=True,
    type=int,
    help="Minumum length of a read after trimming the ends in order not to be discarded",
)
@click.option(
    "-el",
    "--end-length",
    default=150,
    show_default=True,
    type=int,
    help="Window size at either end of the read to be considered as 'ends' for searching",
)
@click.option(
    "--hits/--no-hits",
    "hits_flag",
    default=False,
    is_flag=True,
    type=bool,
    help="""Write the hits identified using the given adapter specifications to TSV. The format is
    standard BLAST outfmt 6 with the following extra columns: read_length (int), discard (bool), trim_l (bool), trim_r (bool)
    """,
)
@click.option(
    "--no-summary",
    default=False,
    type=bool,
    is_flag=True,
    help="Skip writing a summary JSON with the number of hits for each adapter",
)
@click.option(
    "-t",
    "--threads",
    default=1,
    show_default=True,
    type=int,
    help="Number of threads to use for compression",
)
def process_blast(
    blastout: click.Path,
    adapter_yaml: click.Path,
    prefix: str,
    min_length_after_trimming: int,
    end_length: int,
    hits_flag: bool,
    no_summary: bool,
    threads: int,
) -> None:
    """Processes the input blastout file according to the adapter yaml key.

    BLASTOUT: tabular file resulting from a BLAST query of a readset against a BLAST
    database of adapter sequences, run with -outfmt "6 std qlen".

    ADAPTER_YAML: yaml file containing a list with the following fields per adapter:

    \b
    - name: (name of adapter. can be a regular expression)
      discard_middle: True/False (discard read if adapter found in middle)
      discard_end: True/False (discard read if adapter found in end)
      trim_end: True/False (trim read if adapter found in end)
      middle_pident: int (minimum pident requred to identify adapter in middle of read)
      middle_length: int (minimum match length required to identify adapter in middle of read)
      end_pident: int (minimum pident requred to identify adapter in end window)
      end_length: int (minimum match length requred to identify adapter in end window)

    Output:
    By default, writes bgzipped BED to [prefix].bed.gz, and a JSON summary file with raw counts
    of adapter hits detected, counts identified after processing, and the total length of
    removed sequences per adapter to [prefix].summary.json.
    """
    ## Set the prefix. If .blastout is in the filename, we drop
    ## everything including and after it. Otherwise we just drop
    ## the final extension.
    if prefix is None:
        filename = Path(blastout).name
        if ".blastout" in filename:
            prefix = filename.split(".blastout")[0]
        else:
            prefix = Path(blastout).stem

    ## Process BLAST inputs
    results = BlastProcessor(
        blastout, adapter_yaml, end_length, min_length_after_trimming
    )

    ## Write BED file
    logger.info(f"Writing BED file: {prefix}.bed.gz")
    with open(prefix + ".bed.gz", "wb") as f:
        with bgzip.BGZipWriter(f, num_threads=threads) as out:
            out.write(results.bed)

    ## Write output hits if requested
    if hits_flag and results.hits is not None and not results.hits.is_empty():
        logger.info(f"Writing hits file: {prefix}.hits")
        results.generate_hits().write_csv(
            prefix + ".hits", separator="\t", include_header=False
        )

    ## Write summary if not disabled
    if not no_summary:
        logger.info(f"Writing summary file: {prefix}.summary.gz")
        summary = SummariseBlastResults(results)
        summary.generate_summary(prefix + ".summary.json")


@click.command("trim")
@click.argument("bam", type=click.File(mode="rb"))
@click.argument("bed", type=click.Path(exists=True))
@click.argument("outfile", default="-", required=False, type=click.File(mode="wb"))
@click.option(
    "-t",
    "--threads",
    default=1,
    show_default=True,
    type=int,
    help="Number of threads to use for compression",
)
@click.option(
    "--format",
    "-f",
    type=click.Choice(["sam", "bam", "cram"]),
    default="bam",
    show_default=True,
    help="Output file format.",
)
def trim(
    bam: click.File,
    bed: str,
    outfile: click.File,
    threads: int,
    format: str,
) -> None:
    """
    Filter the reads stored in a BAM file using the appropriate BED file produced
    by blastout_to_bed and write to a bgzipped fasta file.

    \b
    BAM: BAM file in which to filter reads
    BED: BED file describing regions of the read set to exclude.
    OUTFILE: Output trimmed BAM file. Defaults to stdout.
    """
    logger.info(f"Filtering {bam.name} using BED file: {bed}")
    logger.info(f"Writing the output to {outfile.name}.")

    filterer = BamTrimmer(
        threads=threads,
        format=format,
    )
    filterer.trim_bam_with_bed(bam, bed, outfile)


cli.add_command(process_blast)
cli.add_command(trim)

if __name__ == "__main__":
    cli()
