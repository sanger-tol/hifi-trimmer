import bgzip
import click
import json
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


@click.group()
def cli():
    """Main entry point for the tool."""
    pass


@click.command("process_blast")
@click.argument("blastout", type=click.File(mode="r"))
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
    "-b",
    "--bam",
    default=None,
    required=False,
    type=click.File(mode="rb"),
    help="If blastout file has no read length field, a BAM file of reads to get read lengths",
)
@click.option(
    "-ml",
    "--min_length_after_trimming",
    default=300,
    type=int,
    help="Minumum length of a read after trimming the ends in order not to be discarded",
)
@click.option(
    "-el",
    "--end_length",
    default=150,
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
    type=int,
    help="Number of threads to use for compression",
)
def process_blast(
    blastout,
    adapter_yaml,
    prefix,
    bam,
    min_length_after_trimming,
    end_length,
    hits_flag,
    no_summary,
    threads,
):
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
    By default, writes BED to [prefix].bed, and a summary file with counts of adapter hits
    detected and filtered to [prefix].summary.
    """
    adapters = read_adapter_yaml(adapter_yaml)
    blast = read_blast(blastout, bam)

    ## Make sure each barcode that matches an adapter matches only one adapter
    check = check_records(blast, adapters)
    if len(check) > 0:
        raise click.ClickException(
            f"Barcodes {', '.join(check)} match to more than one adapter in {adapter_yaml.name}!"
        )

    ## process the blastout file
    hits = match_hits(blast, adapters, end_length)
    actions = determine_actions(hits, end_length, min_length_after_trimming)
    bed = create_bed(actions, end_length)

    if prefix is None:
        prefix = re.sub("\\.blastout.*$", "", blastout.name)

    if hits_flag:
        write_hits(hits).collect().write_csv(
            prefix + ".hits", separator="\t", include_header=False
        )

    if not no_summary:
        with open(prefix + ".summary.json", "w") as f:
            summary = write_summary(blast, hits, actions)
            json.dump(summary, f, indent=4)

    with open(prefix + ".bed.gz", "wb") as f:
        with bgzip.BGZipWriter(f, num_threads=threads) as out:
            out.write(
                bed.collect()
                .write_csv(separator="\t", include_header=False)
                .encode("utf-8")
            )


@click.command("filter_bam")
@click.argument("bam", type=click.File(mode="rb"))
@click.argument("bed", type=click.File(mode="rb"))
@click.argument("outfile", default=None, required=True, type=click.File(mode="wb"))
@click.option(
    "-t",
    "--threads",
    default=1,
    type=int,
    help="Number of threads to use for compression",
)
def filter_bam(bam, bed, outfile, threads):
    """
    Filter the reads stored in a BAM file using the appropriate BED file produced
    by blastout_to_bed and write to a bgzipped fasta file.

    BAM: BAM file in which to filter reads

    BED: BED file describing regions of the read set to exclude.

    OUTFILE: File to write the filtered reads to (bgzipped).
    """
    filters = filter_bam_with_bed(outfile, bam, bed, threads)

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
