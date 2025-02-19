import click
from hifi_trimmer.output import create_bed, filter_bam
from hifi_trimmer.blast import match_hits, determine_actions
from hifi_trimmer.utils import check_records
from hifi_trimmer.read_files import read_adapter_yaml, read_blast


@click.group()
def cli():
    """Main entry point for the tool."""
    pass


@click.command("blastout_to_bed")
@click.argument("blastout", type=click.File(mode="r"))
@click.argument("adapter_yaml", type=click.File(mode="r"))
@click.option(
    "-o",
    "--output",
    default="-",
    required=False,
    type=click.File(mode="wb"),
    help="Output file to write BED to.",
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
    help="Minumum length of a read after trimming the ends in order not to be discarded.",
)
@click.option(
    "-el",
    "--end_length",
    default=150,
    type=int,
    help="Window size at either end of the read to be considered as 'ends' for searching.",
)
@click.option(
    "-hf",
    "--hits_file",
    default=None,
    type=click.File(mode="wb"),
    help="Write the hits identified using the given adapter specifications to the specified TSV.",
)
def blastout_to_bed(
    blastout,
    adapter_yaml,
    output,
    bam,
    min_length_after_trimming,
    end_length,
    hits_file,
):
    """Processes the input blastout file according to the adapter yaml key.

    BLASTOUT: tabular file resulting from BLAST with -outfmt "6 std qlen". If the qlen column
    is missing, lengths can be calculated by passing the --bam option.

    ADAPTER_YAML: yaml file contaning a list with the following fields per adapters:

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
    By default, writes BED to standard output. This can be redirected with the -o/--output option.
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

    if hits_file is not None:
        hits.collect().write_csv(hits_file, separator="\t", include_header=False)

    bed.collect().write_csv(output, separator="\t", include_header=False)


@click.command("filter_bam_to_fasta")
@click.argument("bed", type=click.File(mode="r"))
@click.argument("bam", type=click.File(mode="rb"))
@click.argument("outfile", default=None, required=True, type=click.File(mode="wb"))
def filter_bam_to_fasta(bed, bam, outfile):
    """
    Filter the reads stored in a BAM file using the appropriate BED file produced
    by blastout_to_bed and write to a bgzipped fasta file.

    BED: BED file describing regions of the read set to exclude.
    BAM: BAM file in which to filter reads
    OUTFILE: File to write the filtered reads to (bgzipped).
    """
    filters = filter_bam(outfile, bam, bed)

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


cli.add_command(blastout_to_bed)
cli.add_command(filter_bam_to_fasta)

if __name__ == "__main__":
    cli()
