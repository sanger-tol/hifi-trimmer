from click.testing import CliRunner
from hifi_trimmer.hifi_trimmer import cli
import pathlib
import hashlib
import pytest


def list_test_data(test_dir):
    """
    Get a list of files in the test directory and store them in
    a dict indexed by extension.
    """
    test_data = {}
    for p in pathlib.Path(test_dir).iterdir():
        if p.suffix == ".gz":
            ext = p.suffixes[-2][1:]
        else:
            ext = p.suffix[1:]

        test_data[ext] = str(p.absolute())

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
                testdata["blastout"],
                testdata["yaml"],
            ],
        )

        print(result.stdout)

        md5_bed = md5checksum(pathlib.Path(td) / "test.bed.gz")
        md5_summary = md5checksum(pathlib.Path(td) / "test.summary.json")

    assert result.exit_code == 0
    assert md5_bed == md5checksum(testdata["bed"])
    assert md5_summary == md5checksum(testdata["json"])


@pytest.mark.parametrize("testdata", list_test_dirs())
def test_filter_bam(tmp_path, testdata):
    runner = CliRunner()

    with runner.isolated_filesystem(temp_dir=tmp_path) as td:
        result = runner.invoke(
            cli,
            ["filter_bam", testdata["bam"], testdata["bed"], "test.filtered.fa.gz"],
        )

        md5_fasta = md5checksum(pathlib.Path(td) / "test.filtered.fa.gz")

    assert result.exit_code == 0
    assert md5_fasta == md5checksum(testdata["fa"])


@pytest.mark.parametrize("testdata", list_test_dirs())
def test_filter_bam_fastq(tmp_path, testdata):
    runner = CliRunner()

    with runner.isolated_filesystem(temp_dir=tmp_path) as td:
        result = runner.invoke(
            cli,
            [
                "filter_bam",
                testdata["bam"],
                testdata["bed"],
                "test.filtered.fq.gz",
                "--fastq",
            ],
        )

        md5_fasta = md5checksum(pathlib.Path(td) / "test.filtered.fq.gz")

    assert result.exit_code == 0
    assert md5_fasta == md5checksum(testdata["fq"])
