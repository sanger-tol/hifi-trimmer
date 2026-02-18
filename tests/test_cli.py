import hashlib
import json
import pathlib

import pytest
from click.testing import CliRunner

from hifi_trimmer.hifi_trimmer import cli


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
    """Remove environment-specific paths from summary dict."""
    for key in ["blast_file", "yaml_file"]:
        summary_dict["run_info"].pop(key, None)
    return summary_dict


@pytest.mark.parametrize("testdata", list_test_dirs())
def test_process_blast(tmp_path, testdata):
    runner = CliRunner()

    with runner.isolated_filesystem(temp_dir=tmp_path) as td:
        result = runner.invoke(
            cli,
            [
                "process_blast",
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
def test_filter_bam(tmp_path, testdata):
    runner = CliRunner()

    with runner.isolated_filesystem(temp_dir=tmp_path) as td:
        result = runner.invoke(
            cli,
            [
                "filter_bam",
                testdata["test.bam"],
                testdata["test.bed.gz"],
                "test.filtered.fa.gz",
            ],
        )

        md5_fasta = md5checksum(pathlib.Path(td) / "test.filtered.fa.gz")

    assert result.exit_code == 0
    assert md5_fasta == md5checksum(testdata["test.filtered.fa.gz"])


@pytest.mark.parametrize("testdata", list_test_dirs())
def test_filter_bam_fastq(tmp_path, testdata):
    runner = CliRunner()

    with runner.isolated_filesystem(temp_dir=tmp_path) as td:
        result = runner.invoke(
            cli,
            [
                "filter_bam",
                testdata["test.bam"],
                testdata["test.bed.gz"],
                "test.filtered.fq.gz",
                "--fastq",
            ],
        )

        md5_fasta = md5checksum(pathlib.Path(td) / "test.filtered.fq.gz")

    assert result.exit_code == 0
    assert md5_fasta == md5checksum(testdata["test.filtered.fq.gz"])


@pytest.mark.parametrize("testdata", list_test_dirs())
def test_filter_bam_sam_tags_fastq(tmp_path, testdata):
    runner = CliRunner()

    with runner.isolated_filesystem(temp_dir=tmp_path) as td:
        result = runner.invoke(
            cli,
            [
                "filter_bam",
                testdata["test.bam"],
                testdata["test.bed.gz"],
                "test.filtered.tags.fq.gz",
                "--fastq",
                "--preserve-sam-tags",
            ],
        )

        md5_fasta = md5checksum(pathlib.Path(td) / "test.filtered.tags.fq.gz")

    assert result.exit_code == 0
    assert md5_fasta == md5checksum(testdata["test.filtered.tags.fq.gz"])


@pytest.mark.parametrize("testdata", list_test_dirs())
def test_filter_bam_sam_tags(tmp_path, testdata):
    runner = CliRunner()

    with runner.isolated_filesystem(temp_dir=tmp_path) as td:
        result = runner.invoke(
            cli,
            [
                "filter_bam",
                testdata["test.bam"],
                testdata["test.bed.gz"],
                "test.filtered.tags.fa.gz",
                "--preserve-sam-tags",
            ],
        )

        md5_fasta = md5checksum(pathlib.Path(td) / "test.filtered.tags.fa.gz")

    assert result.exit_code == 0
    assert md5_fasta == md5checksum(testdata["test.filtered.tags.fa.gz"])
