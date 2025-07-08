import bgzip
import click
import importlib.metadata
import json
import polars as pl
import re


from hifi_trimmer.output import (
    create_bed,
    filter_bam_with_bed,
    write_summary,
    write_hits,
)
from hifi_trimmer.blast import match_hits, determine_actions
from hifi_trimmer.utils import check_records
from hifi_trimmer.read_files import read_adapter_yaml, read_blast
from hifi_trimmer.summary import summarise_blast, summarise_hits, summarise_actions


@click.group()
@click.version_option(package_name="hifi_trimmer")
def cli():
    """Main entry point for the tool."""
    pass


@click.command("process_blast")
@click.argument("blastout", type=click.Path(exists=True))
@click.argument("adapter_yaml", type=click.File(mode="r"))
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
    "--min_length_after_trimming",
    default=300,
    show_default=True,
    type=int,
    help="Minumum length of a read after trimming the ends in order not to be discarded",
)
@click.option(
    "-el",
    "--end_length",
    default=150,
    show_default=True,
    type=int,
    help="Window size at either end of the read to be considered as 'ends' for searching",
)
@click.option(
    "-hf",
    "--hits",
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
    help="Skip writing a summary TSV with the number of hits for each adapter",
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
    blastout: str,
    adapter_yaml: click.File,
    prefix: str,
    min_length_after_trimming: int,
    end_length: int,
    hits_flag: bool,
    no_summary: bool,
    threads: int,
) -> None:
    """Processes the input blastout file according to the adapter yaml key.

    BLASTOUT: tabular file resulting from a BLAST query of a readset against a BLAST
    database of adapter sequences, run with -outfmt "6 std qlen". If the qlen column
    is missing, lengths can be calculated by passing the --bam option.

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
    if prefix is None:
        prefix = re.sub("\\.blastout.*$", "", blastout)

    ## Set null outputs
    blast_summary = None
    hits_summary = None
    actions_summary = None
    out_bed = "".encode("utf-8")

    try:
        blast = read_blast(blastout)
        blast_summary = summarise_blast(blast)
        adapters = read_adapter_yaml(adapter_yaml)

        ## Make sure each barcode that matches an adapter matches only one adapter
        check = check_records(blast, adapters)
        if len(check) > 0:
            raise click.ClickException(
                f"Barcodes {', '.join(check)} match to more than one adapter in {adapter_yaml.name}!"
            )

        ## process the blastout file
        hits = match_hits(blast, adapters, end_length).collect()
        hits_summary = summarise_hits(hits)

        ## Check if any hits matched - need to handle empty input!
        if not (hits.is_empty()):
            actions = determine_actions(
                hits.lazy(), end_length, min_length_after_trimming
            ).collect()
            actions_summary = summarise_actions(actions)

            bed = create_bed(actions.lazy(), end_length)

            if hits_flag:
                write_hits(hits).write_csv(
                    prefix + ".hits", separator="\t", include_header=False
                )

            out_bed = (
                bed.collect()
                .write_csv(separator="\t", include_header=False)
                .encode("utf-8")
            )

            click.echo(f"BLAST file {blastout} parsed successfully!")
        else:
            click.echo(
                "WARNING: After filtering, no BLAST hits remain. Writing empty BED output."
            )
    except pl.exceptions.NoDataError:
        click.echo(f"WARNING: {blastout} was empty! Writing empty BED output.")
        blast = None

    if not no_summary:
        with open(prefix + ".summary.json", "w") as f:
            summary = write_summary(blast_summary, hits_summary, actions_summary)
            json.dump(summary, f, indent=4)

    with open(prefix + ".bed.gz", "wb") as f:
        with bgzip.BGZipWriter(f, num_threads=threads) as out:
            out.write(out_bed)


@click.command("filter_bam")
@click.argument("bam", type=click.File(mode="rb"))
@click.argument("bed", type=click.Path(exists=True))
@click.argument("outfile", default=None, required=True, type=click.File(mode="wb"))
@click.option(
    "-f",
    "--fastq",
    is_flag=True,
    default=False,
    help="Write FASTQ instead of FASTA",
)
@click.option(
    "-t",
    "--threads",
    default=1,
    show_default=True,
    type=int,
    help="Number of threads to use for compression",
)
def filter_bam(
    bam: click.File, bed: str, outfile: click.File, threads: int, fastq: bool
) -> None:
    """
    Filter the reads stored in a BAM file using the appropriate BED file produced
    by blastout_to_bed and write to a bgzipped fasta file.

    \b
    BAM: BAM file in which to filter reads
    BED: BED file describing regions of the read set to exclude.
    OUTFILE: File to write the filtered reads to (bgzipped).
    """
    filters = filter_bam_with_bed(outfile, bam, bed, threads, fastq)

    try:
        read = next(filters)
        print(f"Read {read['read']} not processed!")
        for read in filters:
            print(f"Read {read['read']} not processed!")

        raise click.ClickException(
            "Not all entries in the BED file were processed! Did you forget to sort them in the same order as the BAM file?"
        )
    except StopIteration:
        print("Read filtering complete!")


cli.add_command(process_blast)
cli.add_command(filter_bam)

if __name__ == "__main__":
    cli()
