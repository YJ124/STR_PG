from __future__ import annotations

from pathlib import Path

import pytest

from strpg.models import AlleleRecord, LocusRecord
from strpg.registry import AlleleRegistry


def make_locus() -> LocusRecord:
    left = "AACCGGTTAACCGGTTAACCGGTTAACCGG"
    right = "TTGGCCAATTGGCCAATTGGCCAATTGGCC"
    motif = "CAG"
    alleles = []
    for count in (8, 10, 12):
        repeat = motif * count
        alleles.append(
            AlleleRecord(
                f"P1:L{count}", count, left + repeat + right, repeat, motif,
                len(left), len(left) + len(repeat),
                {"GLOBAL": 1 / 3}, {"GLOBAL": 100},
            )
        )
    return LocusRecord(
        "P1", "GRCh38", "chr1", 101, 130, motif, left, right, 12, alleles
    )


@pytest.fixture
def read_selection_fixture(tmp_path: Path):
    locus = make_locus()
    left, right, motif = locus.left_flank, locus.right_flank, locus.motif
    registry = AlleleRegistry({"P1": locus})
    registry_path = tmp_path / "registry.jsonl"
    registry.save(registry_path)

    read1 = left[-12:] + motif * 8 + right[:12]
    read2 = "ACGT" * 12
    quality1 = "I" * len(read1)
    quality2 = "I" * len(read2)
    fq1 = tmp_path / "R1.fastq"
    fq2 = tmp_path / "R2.fastq"
    fq1.write_text(f"@readA/1\n{read1}\n+\n{quality1}\n", encoding="utf-8")
    fq2.write_text(f"@readA/2\n{read2}\n+\n{quality2}\n", encoding="utf-8")
    gaf = tmp_path / "reads.gaf"
    lines = []
    for mate, length in ((1, len(read1)), (2, len(read2))):
        lines.append(
            "\t".join(
                [
                    "readA",
                    str(length),
                    "0",
                    str(length),
                    "+",
                    "P1",
                    "100",
                    "0",
                    str(length),
                    str(length),
                    str(length),
                    "60",
                    "pi:Z:P1",
                    f"mt:i:{mate}",
                ]
            )
        )
    gaf.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {
        "registry": registry,
        "registry_path": registry_path,
        "locus": locus,
        "fq1": fq1,
        "fq2": fq2,
        "gaf": gaf,
        "read1": read1,
        "read2": read2,
    }
