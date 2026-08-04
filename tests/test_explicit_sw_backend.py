from strpg.likelihoods import SWLikelihoodBackend, create_likelihood_backend


def test_explicit_sw_backend():
    assert isinstance(create_likelihood_backend("sw"), SWLikelihoodBackend)
