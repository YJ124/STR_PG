from strpg.likelihoods.base import AlleleLikelihoodBackend
from strpg.read_selection import select_informative_read_records


def test_no_likelihood_called_during_selection(monkeypatch, read_selection_fixture):
    def forbidden(*args, **kwargs):
        raise AssertionError("likelihood backend was called during selection")

    monkeypatch.setattr(AlleleLikelihoodBackend, "score_reads", forbidden)
    fixture = read_selection_fixture
    records = select_informative_read_records(
        sample_id="S1",
        registry=fixture["registry"],
        gaf_path=fixture["gaf"],
        reads1=fixture["fq1"],
        reads2=fixture["fq2"],
    )
    assert len(records) == 2
