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
    md5 = hashlib.md5()
    with open(fname, "rb") as f:
        while chunk := f.read(4096):
            md5.update(chunk)

    return md5.hexdigest()


@pytest.mark.parametrize("testdata", list_test_data())
def test_blastout_to_bed(tmp_path, testdata):
    runner = CliRunner()

    with runner.isolated_filesystem(temp_dir=tmp_path):
        result = runner.invoke(
            cli,
            [
                "blastout_to_bed",
                "--bam",
                testdata["bam"],
                testdata["blastout"],
                testdata["yaml"],
            ],
        )
        result_md5 = hashlib.md5(result.output.encode("utf-8")).hexdigest()
        example_md5 = md5checksum(testdata["bed"])

    assert result.exit_code == 0
    assert result_md5 == example_md5


@pytest.mark.parametrize("testdata", list_test_data())
def test_filter_bam_to_fasta(tmp_path, testdata):
    runner = CliRunner()

    with runner.isolated_filesystem(temp_dir=tmp_path):
        result = runner.invoke(
            cli,
            ["filter_bam_to_fasta", testdata["bed"], testdata["bam"], "output.fa.gz"],
        )
        result_md5 = md5checksum("output.fa.gz")
        example_md5 = md5checksum(testdata["fa"])

    assert result.exit_code == 0
    assert result_md5 == example_md5
