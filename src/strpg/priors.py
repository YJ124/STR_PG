from __future__ import annotations

import math
from itertools import combinations_with_replacement

from .models import AlleleRecord


def parse_population_mix(text: str | None) -> dict[str, float]:
    if not text:
        return {"GLOBAL": 1.0}
    result: dict[str, float] = {}
    for part in text.split(","):
        if not part.strip():
            continue
        pop, value = part.split("=", 1)
        result[pop.strip()] = float(value)
    total = sum(result.values())
    if total <= 0:
        raise ValueError("Population mixture weights must sum to a positive value")
    return {k: v / total for k, v in result.items()}


def mixed_allele_frequencies(
    alleles: list[AlleleRecord],
    population_mix: dict[str, float],
    floor: float = 1e-8,
) -> list[float]:
    values: list[float] = []
    for allele in alleles:
        freq = 0.0
        observed = False
        for pop, weight in population_mix.items():
            if pop in allele.frequencies:
                observed = True
                freq += weight * allele.frequencies[pop]
        if not observed and "GLOBAL" in allele.frequencies:
            freq = allele.frequencies["GLOBAL"]
        values.append(max(floor, freq))
    total = sum(values)
    return [v / total for v in values]


def balding_nichols_prior(pi: float, pj: float, same: bool, theta: float) -> float:
    if not (0.0 <= theta < 1.0):
        raise ValueError("theta must be in [0,1)")
    if same:
        return pi * pi + pi * (1.0 - pi) * theta
    return 2.0 * pi * pj * (1.0 - theta)


def length_similarity(li: int, lj: int, sigma: float) -> float:
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    return math.exp(-((li - lj) ** 2) / (2.0 * sigma * sigma))


def genotype_priors(
    alleles: list[AlleleRecord],
    population_mix: dict[str, float],
    theta: float = 0.01,
    sigma: float = 10.0,
    frequency_floor: float = 1e-8,
    prior_mode: str = "full",
    use_length_smoothing: bool = True,
) -> dict[tuple[int, int], float]:
    """Return normalized priors for unordered diploid genotypes.

    ``prior_mode`` is deliberately explicit for ablation experiments:

    - ``full``: population frequencies with Balding-Nichols correction;
    - ``hwe``: ordinary Hardy-Weinberg equilibrium (theta must be zero);
    - ``uniform``: equal probability for every unordered diploid genotype.

    The default arguments reproduce the production model.
    """
    if prior_mode not in {"full", "hwe", "uniform"}:
        raise ValueError(f"Unsupported prior_mode: {prior_mode}")
    genotypes = list(combinations_with_replacement(range(len(alleles)), 2))
    if not genotypes:
        raise ValueError("At least one allele is required")
    if prior_mode == "uniform":
        value = 1.0 / len(genotypes)
        return {gt: value for gt in genotypes}
    if prior_mode == "hwe" and theta != 0.0:
        raise ValueError("HWE-only prior requires theta=0")
    freqs = mixed_allele_frequencies(alleles, population_mix, floor=frequency_floor)
    raw: dict[tuple[int, int], float] = {}
    effective_theta = 0.0 if prior_mode == "hwe" else theta
    for i, j in genotypes:
        bn = balding_nichols_prior(freqs[i], freqs[j], i == j, effective_theta)
        smooth = (
            length_similarity(alleles[i].repeat_count, alleles[j].repeat_count, sigma)
            if use_length_smoothing
            else 1.0
        )
        raw[(i, j)] = max(1e-300, bn * smooth)
    total = sum(raw.values())
    return {gt: value / total for gt, value in raw.items()}
