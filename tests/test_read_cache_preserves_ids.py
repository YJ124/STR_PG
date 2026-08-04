import csv

from strpg.read_selection import select_informative_read_records, write_read_cache


def test_read_cache_preserves_ids(read_selection_fixture, tmp_path):
    f = read_selection_fixture
    records = select_informative_read_records(
        sample_id="S1", registry=f["registry"], gaf_path=f["gaf"],
        reads1=f["fq1"], reads2=f["fq2"],
    )
    path = write_read_cache(records, tmp_path)[0]
    rows = list(csv.DictReader(path.open(encoding="utf-8"), delimiter="\t"))
    assert [row["read_id"] for row in rows] == ["readA", "readA"]
