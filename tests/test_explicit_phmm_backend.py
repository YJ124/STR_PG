from strpg.likelihoods import PHMMLikelihoodBackend, create_likelihood_backend


def test_explicit_phmm_backend():
    assert isinstance(create_likelihood_backend("phmm"), PHMMLikelihoodBackend)
