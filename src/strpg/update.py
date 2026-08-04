from __future__ import annotations

import csv
import json
from pathlib import Path

from .models import ProvisionalAllele
from .registry import AlleleRegistry
from .utils import open_text


def update_registry(
    registry_in: str | Path,
    genotype_tsv: str | Path,
    registry_out: str | Path,
    population: str,
    min_gq: int = 20,
    base_total: float = 1000.0,
    novel_hypotheses: str | Path | None = None,
    pseudo_count_alpha: float = 0.5,
    sample_id: str | None = None,
    homozygous_only: bool = True,
) -> dict[str, int]:
    registry = AlleleRegistry.load(registry_in)
    stats = {"calls_seen": 0, "calls_used": 0, "allele_copies_added": 0, "novel_accepted": 0, "provisional_stored": 0}

    for locus in registry.loci.values():
        AlleleRegistry.initialize_counts_from_frequencies(locus, base_total=base_total)

    with open_text(genotype_tsv, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            stats["calls_seen"] += 1
            if int(float(row["GQ"])) < min_gq:
                continue
            a1, a2 = int(row["allele1"]), int(row["allele2"])
            if homozygous_only and a1 != a2:
                continue
            locus = registry.loci.get(row["pointer_id"])
            if locus is None:
                continue
            copies = [a1, a2] if not homozygous_only else [a1, a1]
            used = False
            for repeat_count in copies:
                allele = locus.allele_by_count(repeat_count)
                if allele is None:
                    continue
                allele.counts[population] = allele.counts.get(population, 0.0) + 1.0
                stats["allele_copies_added"] += 1
                used = True
            if used:
                stats["calls_used"] += 1

    if novel_hypotheses:
        with open_text(novel_hypotheses, "rt") as handle:
            for raw in handle:
                if not raw.strip():
                    continue
                data = json.loads(raw)
                pointer_id = data.pop("pointer_id")
                hypothesis = ProvisionalAllele.from_dict(data)
                hypothesis.source_sample = sample_id
                if hypothesis.status == "accepted":
                    registry.accept_novel(pointer_id, hypothesis.repeat_count, population, alpha=pseudo_count_alpha, source_sample=sample_id)
                    stats["novel_accepted"] += 1
                else:
                    registry.add_provisional(pointer_id, hypothesis)
                    stats["provisional_stored"] += 1

    for locus in registry.loci.values():
        AlleleRegistry.recompute_frequencies(locus)
    registry.save(registry_out)
    return stats
