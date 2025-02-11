
import click
import polars as pl
import polars.selectors as cs
import pybedtools
import pysam
import re
import sys
import yaml

def read_adapter_yaml(yaml_path: str) -> pl.DataFrame:
    with open(yaml_path) as stream:
        try:
            adapters = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    return pl.DataFrame(adapters).lazy()

def read_lengths(filtered_reads: set, bam: str) -> pl.DataFrame:
    """Count the lengths of reads in set filtered_reads using bamfile bam"""
    lengths = []
    with pysam.AlignmentFile(bam, "rb", check_sq=False, require_index=False) as b:
        for read in b.fetch(until_eof=True):
            if read.query_name in filtered_reads:
                lengths.append({
                    'qseqid': read.query_name,
                    'read_length': read.query_length
                })

    return pl.DataFrame(lengths)

def match_hits(blastout: pl.DataFrame, adapter_df: str, end_length: int) -> pl.DataFrame:
    """Determine whether each blast hit matches sufficient criteria for discard/trimming
    
    Keyword arguments:
    blastout: pl.DataFrame containing the results of a blast search with read lengths in column 13
    adapter_key: pl.DataFrame of adapters and their filtering/trimming parameters
    end_length: window at either end of read which is examined for trimming
    """
    discard = (
        ((pl.col("discard_middle")).and_(
            (pl.col("qstart") > end_length),
            (pl.col("qend") < pl.col("read_length") - end_length),
            (pl.col("pident") >= pl.col("middle_pident")),
            (pl.col("length") >= pl.col("middle_length"))
        )) |
        ((pl.col("discard_end")).and_(
            (pl.col("qend") < end_length),
            (pl.col("qstart") > pl.col("read_length") - end_length),
            (pl.col("pident") >= pl.col("end_pident")),
            (pl.col("length") >= pl.col("end_length"))
        ))
    )

    trim_l = (
        (pl.col("trim_end"))
        .and_(
            pl.col("qend") <= end_length,
            (pl.col("pident") >= pl.col("end_pident")),
            (pl.col("length") >= pl.col("end_length"))
        )
    )

    trim_r = (
        (pl.col("trim_end"))
        .and_(
            pl.col("qstart") >= pl.col("read_length") - end_length,
            (pl.col("pident") >= pl.col("end_pident")),
            (pl.col("length") >= pl.col("end_length"))
        )
    )
    
    return (blastout
        .join(adapter_df, how="cross", on=None, maintain_order="left")
        .filter(
            pl.struct(["sseqid", "barcode"]).map_elements(
                lambda s: bool(re.search(s["barcode"], s["sseqid"])), return_dtype=pl.Boolean
            )
        )
        .with_columns(
            discard=discard,
            trim_l=trim_l,
            trim_r=trim_r
        )
        .filter(pl.any_horizontal("discard", "trim_l", "trim_r"))

    )

def create_bed(blastout: pl.DataFrame, end_length: int, min_length: int) -> pl.DataFrame:
    blastout = (blastout
        .with_columns(
            (pl.col("read_length") - ((pl.col("trim_l").any() + pl.col("trim_r").any()) * end_length)).over("qseqid").alias("read_length_after_trimming"),
        )
        .unpivot(
            cs.starts_with("trim_"),
            index=~cs.starts_with("trim_"),
            variable_name="trim_where",
            value_name="trim"
        )
    )

    discard = (blastout
        .filter((pl.col("discard").any() | pl.col("read_length_after_trimming").lt(min_length).any()).over("qseqid"))
        .group_by("qseqid")
        .agg(
            start = pl.lit(0, dtype=pl.Int64),
            end = pl.col("read_length").first(),
            reason = pl.concat_str(pl.lit("discard:"), pl.col("sseqid").first())
        )
    )

    trim = (
        blastout
        .filter(~(pl.col("discard").any() | pl.col("read_length_after_trimming").lt(min_length).any()).over("qseqid"))
        .filter(pl.col("trim"))
        .with_columns(
            pl.when(
                (pl.col("trim_where") == "trim_l")
            ).then(pl.struct(
                pl.lit(0).alias("start"),
                pl.lit(end_length).alias("end")
                )
            ).otherwise(pl.struct(
                (pl.col("read_length") - end_length).alias("start"),
                pl.col("read_length").alias("end")
                )
            ).struct.field(["start", "end"])
        )
        .group_by("qseqid", "trim_where")
        .agg(
            pl.col("start").min().alias("start"),
            pl.col("end").max().alias("end"),
            pl.concat_str(
                pl.lit("trim:"),
                pl.col("sseqid").get(pl.col("qend").arg_max())
            ).alias("reason")
        )
        .select(["qseqid", "start", "end", "reason"])
    )

    result = pl.concat([discard, trim])
    return result

@click.group()
def cli():
    """Main entry point for the tool."""
    pass

@click.command("blastout_to_bed")
@click.argument('blastout', type=click.Path(exists=True))
@click.argument('adapter_yaml', type=click.Path(exists=True))
@click.argument('outfile', default=None, required=False, type=click.File("w"))
@click.option('--bam', default=None, required=False, type=click.Path(exists=True), 
    help="If blastout file has no read length field, a BAM file of reads to get read lengths")
@click.option("--min_length_after_trimming", default=300, type=int,
    help="Minumum length of a read after trimming the ends in order not to be discarded.")
@click.option("--end_length", default=150, type=int,
    help = "Window size at either end of the read to be considered as 'ends' for searching.")
def blastout_to_bed(blastout, adapter_yaml, bam, outfile, min_length_after_trimming, end_length):
    """Processes the input blastout file according to the adapter yaml key.
    
    BLASTOUT: tabular file resulting from BLAST with -outfmt "6 std qlen". If the qlen column
    is missing, lengths can be calculated by passing the --bam option.
    ADAPTER_YAML: yaml file contaning a list with the following fields per barcode:
        - name: name of adapter. can be a regular expression
          discard_middle: True/False - discard read if adapter found in middle
          discard_end: True/False  - discard read if adapter found in end
          trim_end: True/False - trim read if adapter found in end
          middle_pident: int - minimum pident requred to identify adapter in middle of read
          middle_length: int - minimum match length required to identify adapter in middle of read
          end_pident: int - minimum pident requred to identify adapter in end window
          end_length: int - minimum match length requred to identify adapter in end window
    OUTFILE: (optional) - file to write results to - prints to STDOUT if not supplied
    """ 
    ## read in data
    adapters = read_adapter_yaml(adapter_yaml)
    blast = pl.scan_csv(blastout, has_header=False, separator='\t')

    ## rename blast columns based on number of input columns
    blast_ncol = len(blast.head(1).collect().columns)

    column_names = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]
    if blast_ncol == 13:
        column_names = column_names + ["read_length"]

    blast = (blast
        .rename({old: new for old, new in zip(blast.collect_schema().names(), column_names)})
        ## reorder start-end to correct strandedness as we're not interested here
        .with_columns(
            pl.when(
                pl.col("qstart") > pl.col("qend")
            )
            .then(pl.struct(
                start = "qend",
                end = "qstart",
            ))
            .otherwise(pl.struct(["qstart", "qend"]))
            .struct.field(["qstart", "qend"])
        )
    )
    
    ## if read lengths not provided but bam is - get them from the bam
    if(blast_ncol == 12 and bam is not None):
        filtered_reads = set(blast.select("qseqid").unique().collect()["qseqid"].to_list())
        lengths = read_lengths(filtered_reads, bam)
        blast = blast.join(lengths.lazy(), on = "qseqid", how = "left", suffix="")
    else:
        sys.exit("Error: BLAST file only has 12 columns, but no BAM file has been provided to count!")
    
    ## process the blastout file
    blast = match_hits(blast, adapters, end_length)
    result = create_bed(blast, end_length, min_length_after_trimming).collect()

    if outfile:
        click.echo(f"writing bed file to {outfile}")
        result.write_csv(file = outfile, separator="\t", include_header=False)
    else:
        click.echo(result.write_csv(separator="\t", include_header=False))

    return 0


@click.command("filter_bam_to_fasta")
@click.argument('bed', type=click.Path(exists=True))
@click.argument('bam', type=click.Path(exists=True))
@click.argument('outfile', default=None, required=False, type=click.File("w"))
def filter_bam_to_fasta(bed, bam, outfile):
    """
    Filter the reads stored in a BAM file using the appropriate BED file produced 
    by blastout_to_bed.

    BED: BED file describing regions of the read set to exclude.
    BAM: BAM file in which to filter reads
    OUTFILE: The output fasta file
    """
    ## Create genome file - needed for complement
    lengths = []
    with pysam.AlignmentFile(bam, "rb", check_sq=False, require_index=False) as bam: 
        for read in bam.fetch(until_eof=True):
            lengths.append((read.name, read.length))

    genome = pybedtools.genome.Genome(dict(lengths))

    ## Load BED file and complement it
    bedtool = pybedtools.BedTool(bed).complement(g = genome)

    

cli.add_command(blastout_to_bed)
cli.add_command(filter_bam_to_fasta)

if __name__ == '__main__':
    cli()
