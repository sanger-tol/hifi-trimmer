from tempfile import TemporaryDirectory

import duckdb

from hifi_trimmer import read_files


class HiFiTrimmerDB:
    def __init__(self, dbfile=None):
        self.tmp_dir = TemporaryDirectory()
        if not dbfile:
            dbfile = f"{self.tmp_dir.name}/hifi_trimmer.duckdb"
        self.conn = duckdb.connect(dbfile)
        self._setup_duckdb()

    def _setup_duckdb(self):
        """
        Farm friendly settings. Limit memory usage, and use /tmp for temporary
        storage, which should be fast SSD. May also want to limit threads here.
        """
        self.conn.execute("SET memory_limit = ?", ("4GB",))
        self.conn.execute("SET temp_directory = ?", (str(self.tmp_dir),))

    def read_blast(self, blast_path: str, table_name="blastout"):
        self.conn.execute(
            f"""
            CREATE TABLE {table_name} AS
            SELECT
                qseqid
              , sseqid
              , pident
              , length
              , mismatch
              , gapopen
              , if(qstart > qend, qend, qstart) AS qstart
              , if(qstart > qend, qstart, qend) AS qend
              , sstart
              , send
              , evalue
              , bitscore
              , read_length
            FROM read_csv(?, columns = {{
                  'qseqid': 'VARCHAR'
                , 'sseqid': 'VARCHAR'
                , 'pident': 'FLOAT'
                , 'length': 'INTEGER'
                , 'mismatch': 'INTEGER'
                , 'gapopen': 'INTEGER'
                , 'qstart': 'INTEGER'
                , 'qend': 'INTEGER'
                , 'sstart': 'INTEGER'
                , 'send': 'INTEGER'
                , 'evalue': 'FLOAT'
                , 'bitscore': 'FLOAT'
                , 'read_length': 'INTEGER'
              }}
            )
            """, (blast_path,)
        )

    def read_adapter_yaml(self, yaml_file: str, table_name="adapters"):
        """
        Creates a table from the Polars dataframe
        """
        aptrs = read_files.read_adapter_yaml(yaml_file)  # noqa: F841
        print(f"{aptrs = }")
        self.conn.execute(f"CREATE TABLE {table_name} AS FROM aptrs")

