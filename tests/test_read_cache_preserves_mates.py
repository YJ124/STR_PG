from strpg.read_selection import select_informative_read_records


def test_read_cache_preserves_mates(read_selection_fixture):
    f = read_selection_fixture
    records = select_informative_read_records(
        sample_id="S1", registry=f["registry"], gaf_path=f["gaf"],
        reads1=f["fq1"], reads2=f["fq2"],
    )
    assert [record.mate for record in records] == [1, 2]
