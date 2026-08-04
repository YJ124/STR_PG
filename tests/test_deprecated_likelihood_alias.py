import pytest

from strpg.cli import main


def test_deprecated_alias_emits_warning_before_io():
    with pytest.warns(FutureWarning, match="deprecated"):
        with pytest.raises(FileNotFoundError):
            main(
                [
                    "genotype", "--registry", "missing", "--gaf", "missing",
                    "--fq1", "missing", "--out", "missing",
                    "--likelihood-model", "alignment_softmax",
                ]
            )
