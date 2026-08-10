import hashlib
import json
import pathlib

import pysam
import pytest
from click.testing import CliRunner

from hifi_trimmer.cli import cli


def list_test_data(test_dir):
    """
    Get a list of files in the test directory and store them in
    a dict indexed by full filename.
    """
    test_data = {}
    for p in pathlib.Path(test_dir).iterdir():
        # Use the full filename as the key
        key = p.name
        test_data[key] = str(p.absolute())

    return test_data


def bam_to_records(path):
    """Extract (name, sequence, qualities) tuples from a BAM/SAM/CRAM file."""
    with pysam.AlignmentFile(str(path), check_sq=False) as f:
        return [
            (
                r.query_name,
                r.query_sequence,
                r.query_qualities.tobytes() if r.query_qualities else None,
            )
            for r in f.fetch(until_eof=True)
        ]


def fastx_to_records(path):
    """Extract (name, sequence, qualities) tuples from a FASTA/FASTQ file.
    
    Returns None for qualities in FASTA files.
    Qualities are returned as a string for FASTQ files.
    """
    records = []
    with pysam.FastxFile(str(path)) as f:
        for r in f:
            quality = r.quality if r.quality else None
            records.append((r.name, r.sequence, quality))
    return records


def list_test_dirs():
    """
    List all the data dirs in the test data directory.
    """
    data_dir = pathlib.Path(__file__).parent / "data"

    data = []
    for p in pathlib.Path(data_dir).iterdir():
        if p.is_dir():
            data.append(list_test_data(p))

    return data


def md5checksum(fname):
    """
    Get the MD5 checksum for a file.
    """
    md5 = hashlib.md5()
    with open(fname, "rb") as f:
        while chunk := f.read(4096):
            md5.update(chunk)

    return md5.hexdigest()


def clean_summary_for_comparison(summary_dict):
    """Remove environment-specific paths and version info from summary dict."""
    for key in ["blast_file", "yaml_file", "hifi_trimmer"]:
        summary_dict["run_info"].pop(key, None)
    return summary_dict


@pytest.mark.parametrize("testdata", list_test_dirs())
def test_process_blast(tmp_path, testdata):
    runner = CliRunner()

    with runner.isolated_filesystem(temp_dir=tmp_path) as td:
        result = runner.invoke(
            cli,
            [
                "process-blast",
                "--prefix",
                "test",
                testdata["test.blastout.gz"],
                testdata["test.yaml"],
            ],
        )

        print(result.stdout)

        md5_bed = md5checksum(pathlib.Path(td) / "test.bed.gz")
        summary_dict = clean_summary_for_comparison(
            json.loads((pathlib.Path(td) / "test.summary.json").read_text())
        )

    comparison_summary = clean_summary_for_comparison(
        json.loads(pathlib.Path(testdata["test.summary.json"]).read_text())
    )

    assert result.exit_code == 0
    assert md5_bed == md5checksum(testdata["test.bed.gz"])
    assert summary_dict == comparison_summary


@pytest.mark.parametrize("testdata", list_test_dirs())
@pytest.mark.parametrize("fmt", ["sam", "bam", "cram"])
def test_trim_fmt(tmp_path, testdata, fmt):
    runner = CliRunner()
    outfile = f"test.filtered.{fmt}"
    with runner.isolated_filesystem(temp_dir=tmp_path) as td:
        result = runner.invoke(
            cli,
            [
                "trim",
                "-f",
                fmt,
                "-o",
                outfile,
                testdata["test.bam"],
                testdata["test.bed.gz"],
            ],
        )
        assert result.exit_code == 0
        assert bam_to_records(pathlib.Path(td) / outfile) == bam_to_records(
            testdata[outfile]
        )


@pytest.mark.parametrize("testdata", list_test_dirs())
def test_trim_fasta_from_fasta(tmp_path, testdata):
    """Test trimming FASTA input and writing FASTA output."""
    runner = CliRunner()
    outfile = "test.filtered.fa.gz"
    with runner.isolated_filesystem(temp_dir=tmp_path) as td:
        result = runner.invoke(
            cli,
            [
                "trim",
                "-f",
                "fasta",
                "-o",
                outfile,
                testdata["test.fasta.gz"],
                testdata["test.bed.gz"],
            ],
        )
        assert result.exit_code == 0
        assert fastx_to_records(pathlib.Path(td) / outfile) == fastx_to_records(
            testdata[outfile]
        )


@pytest.mark.parametrize("testdata", list_test_dirs())
def test_trim_fasta_from_fastq(tmp_path, testdata):
    """Test trimming FASTQ input and writing FASTA output."""
    runner = CliRunner()
    outfile = "test.filtered.fa.gz"
    with runner.isolated_filesystem(temp_dir=tmp_path) as td:
        result = runner.invoke(
            cli,
            [
                "trim",
                "-f",
                "fasta",
                "-o",
                outfile,
                testdata["test.fastq.gz"],
                testdata["test.bed.gz"],
            ],
        )
        assert result.exit_code == 0
        assert fastx_to_records(pathlib.Path(td) / outfile) == fastx_to_records(
            testdata[outfile]
        )


@pytest.mark.parametrize("testdata", list_test_dirs())
def test_trim_fastq_from_fastq(tmp_path, testdata):
    """Test trimming FASTQ input and writing FASTQ output."""
    runner = CliRunner()
    outfile = "test.filtered.fq.gz"
    with runner.isolated_filesystem(temp_dir=tmp_path) as td:
        result = runner.invoke(
            cli,
            [
                "trim",
                "-f",
                "fastq",
                "-o",
                outfile,
                testdata["test.fastq.gz"],
                testdata["test.bed.gz"],
            ],
        )
        assert result.exit_code == 0
        assert fastx_to_records(pathlib.Path(td) / outfile) == fastx_to_records(
            testdata[outfile]
        )


@pytest.mark.parametrize("testdata", list_test_dirs())
def test_trim_fasta_from_bam(tmp_path, testdata):
    """Test trimming BAM input and writing FASTA output."""
    runner = CliRunner()
    outfile = "test.filtered.fa.gz"
    with runner.isolated_filesystem(temp_dir=tmp_path) as td:
        result = runner.invoke(
            cli,
            [
                "trim",
                "-f",
                "fasta",
                "-o",
                outfile,
                testdata["test.bam"],
                testdata["test.bed.gz"],
            ],
        )
        assert result.exit_code == 0
        assert fastx_to_records(pathlib.Path(td) / outfile) == fastx_to_records(
            testdata[outfile]
        )


@pytest.mark.parametrize("testdata", list_test_dirs())
def test_trim_fastq_from_bam(tmp_path, testdata):
    """Test trimming BAM input and writing FASTQ output."""
    runner = CliRunner()
    outfile = "test.filtered.fq.gz"
    with runner.isolated_filesystem(temp_dir=tmp_path) as td:
        result = runner.invoke(
            cli,
            [
                "trim",
                "-f",
                "fastq",
                "-o",
                outfile,
                testdata["test.bam"],
                testdata["test.bed.gz"],
            ],
        )
        assert result.exit_code == 0
        assert fastx_to_records(pathlib.Path(td) / outfile) == fastx_to_records(
            testdata[outfile]
        )
