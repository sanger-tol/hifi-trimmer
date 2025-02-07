
import click
import csv
import polars as pl
import pysam
import re

def tsv_to_dict(filename):
    with open(filename, mode='r', newline='', encoding='utf-8') as file:
        reader = csv.DictReader(file, delimiter='\t')
        data = []
        for row in reader:
            for key, value in row.items():
                if value == '':
                    row[key] = None
            data.append(row)
    return data

@click.group()
def cli():
    """Main entry point for the tool."""
    pass

@click.command("process_blast")
@click.argument('blastout', type=click.file('r'))
@click.argument('adapter_key', type=click.file('r'))
@click.argument('bam', default=None, required=False)
def process_blast(blastout, adapter_key, lengths):
    """Processes the input blastout file according to the adapter key"""
    key_dict = tsv_to_dict(adapter_key)
    blast = pl.scan_csv(blastout, 
        has_header=False
    )

    if(blast.shape[1] == 12 and bam is not None):
        filtered_reads = set(blast.select("column_1").unique().to_list())

        lengths = []
        with pysam.AlignmentFile(bam, "rb") as bam:
            for read in bam.fetch():
                if read.query_name in filtered_reads:
                    results.append({
                        'column_1': read.query_name,
                        'length': read.query_length
                    })

        lengths = pl.DataFrame(lengths)
        blast = blast.join(lengths, how = "left")
    
    print(blast)

# @click.command()
# @click.argument('num1', type=int)
# @click.argument('num2', type=int)
# def add_numbers(num1, num2):
#     """Adds two numbers."""
#     result = num1 + num2
#     click.echo(f"The sum is: {result}")

# Register subcommands with the main CLI group
cli.add_command(process_blast)
# cli.add_command(add_numbers)

if __name__ == '__main__':
    cli()
