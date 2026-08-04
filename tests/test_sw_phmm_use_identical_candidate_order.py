from strpg.likelihoods import create_likelihood_backend


def test_sw_phmm_use_identical_candidate_order(read_selection_fixture):
    f = read_selection_fixture
    candidates = f["locus"].alleles
    reads = [f["read1"], f["read2"]]
    sw = create_likelihood_backend("sw").score_reads(reads, candidates, f["locus"])
    phmm = create_likelihood_backend("phmm").score_reads(
        reads, candidates, f["locus"]
    )
    assert sw.shape == phmm.shape
    assert [a.repeat_count for a in candidates] == [8, 10, 12]
