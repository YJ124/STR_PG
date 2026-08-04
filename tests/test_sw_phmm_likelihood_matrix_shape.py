from strpg.likelihoods import create_likelihood_backend


def test_sw_phmm_likelihood_matrix_shape(read_selection_fixture):
    f = read_selection_fixture
    reads = [f["read1"], f["read2"]]
    alleles = f["locus"].alleles
    for name in ("sw", "phmm"):
        assert create_likelihood_backend(name).score_reads(
            reads, alleles, f["locus"]
        ).shape == (2, 3)
