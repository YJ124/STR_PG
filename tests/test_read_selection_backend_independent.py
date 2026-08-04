from strpg.read_selection import select_informative_read_records


def test_read_selection_backend_independent(read_selection_fixture):
    fixture = read_selection_fixture
    records = select_informative_read_records(
        sample_id="S1",
        registry=fixture["registry"],
        gaf_path=fixture["gaf"],
        reads1=fixture["fq1"],
        reads2=fixture["fq2"],
    )
    assert [record.read_id for record in records] == ["readA", "readA"]
    assert all(record.selected for record in records)
