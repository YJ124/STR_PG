from strpg.cli import build_parser
from strpg.defaults import DEFAULT_LIKELIHOOD_BACKEND


def test_default_backend_is_sw():
    args = build_parser().parse_args(
        [
            "genotype",
            "--registry", "registry.jsonl",
            "--gaf", "reads.gaf",
            "--fq1", "R1.fastq",
            "--out", "calls.tsv",
        ]
    )
    assert DEFAULT_LIKELIHOOD_BACKEND == "sw"
    assert args.likelihood_backend == "sw"
