
import click
import polars as pl
import polars.selectors as cs
import pysam
import re
import sys

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

def determine_read_actions(blastout: pl.DataFrame, adapter_key: str, end_length: int) -> pl.DataFrame:
    """Determine whether to discard or trim each read based on adapter position
    
    Keyword arguments:
    blastout: pl.DataFrame containing the results of a blast search with read lengths in column 13
    adapter_key: path to TSV file of adapters and their filtering/trimming parameters
    end_length: window at either end of read which is examined for trimming
    """
    barcodes = pl.scan_csv(adapter_key, separator="\t")

    discard = (
        ((pl.col("dm") == 1).and_(
            (pl.col("pident") >= pl.col("mp")),
            (pl.col("length") >= pl.col("ml"))
        )) |
        ((pl.col("de") == 1).and_(
            (pl.col("pident") >= pl.col("ep")),
            (pl.col("length") >= pl.col("el"))
        ))
    )

    trim_l = (
        (pl.col("te") == 1)
        .and_(
            pl.col("qend") <= end_length,
            (pl.col("pident") >= pl.col("ep")),
            (pl.col("length") >= pl.col("el"))
        )
    )

    trim_r = (
        (pl.col("te") == 1)
        .and_(
            pl.col("qstart") >= pl.col("read_length") - end_length,
            (pl.col("pident") >= pl.col("ep")),
            (pl.col("length") >= pl.col("el"))
        )
    )
    
    return (blastout
        .join(barcodes, how="cross", on=None)
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
    )

def discard_reads(blastout: pl.DataFrame, reason: pl.Expr) -> pl.DataFrame:
    """Takes a dataframe of blast entries of reads to discard and a pl.Expr expressing
    the reason to discard. Returns a pl.DataFrame in BED format for output.
    """
    return (blastout
        .group_by(pl.col("qseqid"))
        .agg(
            pl.lit(0, dtype=pl.Int64()).alias("start"),
            pl.col("read_length").first().alias("end"),
            pl.concat_str(
                pl.lit("discard:"),
                reason
            ).alias("reason")
        )
    )

def trim_reads(blastout: pl.DataFrame, end_length: int, min_length: int) -> pl.DataFrame:
    """Takes a dataframe of blast entries of reads to trim, the length of the window at either
    end of the reads to find hits in, and the minimum length a read must have after trimming 
    to be kept.
    
    Returns a pl.DataFrame in BED format for output.
    """

    blastout_tooshort = discard_reads(
        blastout.filter(pl.col("read_length") - ((pl.col("trim_l") + pl.col("trim_r")) * end_length) < min_length),
        pl.lit("trimmed_tooshort")
    )

    blastout_trimmed = (blastout
        .filter(pl.col("read_length") - ((pl.col("trim_l") + pl.col("trim_r")) * end_length) >= min_length)
        .unpivot(
            cs.starts_with("trim_"),
            index=~cs.starts_with("trim_"),
            variable_name="trim_where",
            value_name="trim"
        )
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
            ).struct.field(["start", "end"]),
            pl.lit("trim").alias("reason")
        )
        .group_by(pl.col("qseqid"), pl.col("trim_where"))
        .agg(
            pl.col("start").min().alias("start"),
            pl.col("end").max().alias("end"),
            pl.concat_str(
                pl.lit("trim:"),
                pl.col("sseqid").unique().str.join(",")
            ).alias("reason")
        )
        .select(["qseqid", "start", "end", "reason"])
    )

    return pl.concat([blastout_tooshort, blastout_trimmed])

@click.group()
def cli():
    """Main entry point for the tool."""
    pass

@click.command("blastout_to_bed")
@click.argument('blastout', type=click.Path(exists=True))
@click.argument('adapter_key', type=click.Path(exists=True))
@click.argument('outfile', default=None, required=False, type=click.File("w"))
@click.option('--bam', default=None, required=False, type=click.Path(exists=True), 
    help="If blastout file has no read length field, a BAM file of reads to get read lengths")
@click.option("--min_length_after_trimming", default=300, type=int,
    help="Minumum length of a read after trimming the ends in order not to be discarded.")
@click.option("--end_length", default=150, type=int,
    help = "Window size at either end of the read to be considered as 'ends' for searching.")
def blastout_to_bed(blastout, adapter_key, bam, outfile, min_length_after_trimming, end_length):
    """Processes the input blastout file according to the adapter key""" 
    blast = pl.scan_csv(blastout, 
        has_header=False,
        separator='\t', 
    )
    blast_ncol = len(blast.head(1).collect().columns)

    column_names = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]
    if blast_ncol == 13:
        column_names = column_names + ["read_length"]

    ## Name columns and sort start and end of query
    blast = (blast
        .rename({old: new for old, new in zip(blast.collect_schema().names(), column_names)})
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
    
    if(blast_ncol == 12 and bam is not None):
        filtered_reads = set(blast.select("qseqid").unique().collect()["qseqid"].to_list())
        lengths = read_lengths(filtered_reads, bam)
        blast = blast.join(lengths.lazy(), on = "qseqid", how = "left", suffix="")
    else:
        sys.exit("Error: BLAST file only has 12 columns, but no BAM file has been provided to count!")
    
    blast = determine_read_actions(blast, adapter_key, end_length)
    discarded_reads = discard_reads(
        blast.filter(pl.col("discard")),
        pl.col("sseqid").unique().str.join(",")
    )
    trimmed_reads = trim_reads(
        blast.filter(pl.col("discard").eq(False).and_((pl.col("trim_l")).or_(pl.col("trim_r")))),
        end_length, 
        min_length_after_trimming
    )

    result = pl.concat([discarded_reads, trimmed_reads]).collect()

    if outfile:
        click.echo(f"writing bed file to {outfile}")
        result.write_csv(file = outfile, separator="\t", include_header=False)
    else:
        click.echo(result.write_csv(separator="\t", include_header=False))

    return 0

cli.add_command(blastout_to_bed)

if __name__ == '__main__':
    cli()
