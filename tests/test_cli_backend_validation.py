import pytest

from strpg.cli import build_parser, main


def test_cli_rejects_unknown_backend():
    with pytest.raises(SystemExit):
        build_parser().parse_args(
            [
                "genotype", "--registry", "r", "--gaf", "g", "--fq1", "f",
                "--out", "o", "--likelihood-backend", "unknown",
            ]
        )


def test_cli_rejects_conflicting_backend_options():
    with pytest.warns(FutureWarning):
        with pytest.raises(SystemExit, match="Conflicting"):
            main(
                [
                    "genotype", "--registry", "r", "--gaf", "g", "--fq1", "f",
                    "--out", "o", "--likelihood-backend", "phmm",
                    "--likelihood-model", "sw",
                ]
            )
