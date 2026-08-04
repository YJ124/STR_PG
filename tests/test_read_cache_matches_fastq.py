import hashlib

from strpg.read_selection import select_informative_read_records


def test_read_cache_matches_fastq(read_selection_fixture):
    f = read_selection_fixture
    records = select_informative_read_records(
        sample_id="S1", registry=f["registry"], gaf_path=f["gaf"],
        reads1=f["fq1"], reads2=f["fq2"],
    )
    assert records[0].sequence == f["read1"]
    assert records[1].sequence == f["read2"]
    assert records[0].sequence_sha256 == hashlib.sha256(
        f["read1"].encode("ascii")
    ).hexdigest()
