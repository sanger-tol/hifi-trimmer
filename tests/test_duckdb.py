import pathlib

from hifi_trimmer.duckdb import HiFiTrimmerDB


def test_duckdb_creation():
    hits_dir = pathlib.Path(__file__).parent / "data" / "hits"
    hits_blast = hits_dir / "test.blastout.gz"
    hits_yaml = hits_dir / "test.yaml"
    db = HiFiTrimmerDB()
    db.read_blast(str(hits_blast))
    assert db.conn.execute("SELECT COUNT(*) FROM blastout").fetchone() == (363,)
    db.read_adapter_yaml(hits_yaml.open())
    assert db.conn.execute("SELECT COUNT(*) FROM adapters").fetchone() == (2,)
