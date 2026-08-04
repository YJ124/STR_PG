from strpg.read_selection import select_informative_read_records


def test_sw_phmm_use_identical_read_ids(read_selection_fixture):
    f = read_selection_fixture
    records = select_informative_read_records(
        sample_id="S1", registry=f["registry"], gaf_path=f["gaf"],
        reads1=f["fq1"], reads2=f["fq2"],
    )
    ids_for_sw = [f"{r.read_id}/{r.mate}" for r in records if r.selected]
    ids_for_phmm = [f"{r.read_id}/{r.mate}" for r in records if r.selected]
    assert ids_for_sw == ids_for_phmm
