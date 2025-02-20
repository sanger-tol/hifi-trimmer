from click.testing import CliRunner
from hifi_trimmer.hifi_trimmer import cli
import pathlib
import hashlib
import pytest


def list_test_data():
    """
    Get a list of files in the test directory and store them in
    a dict indexed by extension.
    """
    data_dir = pathlib.Path(__file__).parent / "data"

    test_data = {}
    for p in pathlib.Path(data_dir).iterdir():
        if p.suffix == ".gz":
            ext = p.suffixes[-2][1:]
        else:
            ext = p.suffix[1:]

        test_data[ext] = str(p.absolute())

    return [test_data]


def md5checksum(fname):
    """
    Get the MD5 checksum for a file.
    """
    md5 = hashlib.md5()
    with open(fname, "rb") as f:
        while chunk := f.read(4096):
            md5.update(chunk)

    return md5.hexdigest()


@pytest.mark.parametrize("testdata", list_test_data())
def test_process_blast(tmp_path, testdata):
    runner = CliRunner()

    with runner.isolated_filesystem(temp_dir=tmp_path) as td:
        result = runner.invoke(
            cli,
            [
                "process_blast",
                "--prefix",
                "test",
                "--hits",
                "--bam",
                testdata["bam"],
                testdata["blastout"],
                testdata["yaml"],
            ],
        )

        md5_bed = md5checksum(pathlib.Path(td) / "test.bed.gz")
        md5_summary = md5checksum(pathlib.Path(td) / "test.summary.json")
        md5_hits = md5checksum(pathlib.Path(td) / "test.hits")

    assert result.exit_code == 0
    assert md5_bed == md5checksum(testdata["bed"])
    assert md5_summary == md5checksum(testdata["json"])
    assert md5_hits == md5checksum(testdata["hits"])


@pytest.mark.parametrize("testdata", list_test_data())
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
