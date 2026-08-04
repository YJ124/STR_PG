from strpg.read_selection import select_informative_read_records, write_read_cache


def test_read_cache_reproducible(read_selection_fixture, tmp_path):
    f = read_selection_fixture
    records = select_informative_read_records(
        sample_id="S1", registry=f["registry"], gaf_path=f["gaf"],
        reads1=f["fq1"], reads2=f["fq2"],
    )
    first = write_read_cache(records, tmp_path / "a")[0].read_bytes()
    second = write_read_cache(records, tmp_path / "b")[0].read_bytes()
    assert first == second
