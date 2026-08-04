from __future__ import annotations

import random
from pathlib import Path

from .registry import AlleleRegistry
from .utils import open_text, reverse_complement


def _mutate(seq: str, error_rate: float, rng: random.Random) -> str:
    bases = "ACGT"
    out = []
    for base in seq:
        if rng.random() < error_rate and base in bases:
            out.append(rng.choice([b for b in bases if b != base]))
        else:
            out.append(base)
    return "".join(out)


def simulate_reads(
    registry_path: str | Path,
    pointer_id: str,
    weights_text: str,
    out1: str | Path,
    out2: str | Path | None = None,
    coverage: int = 30,
    read_length: int = 100,
    fragment_mean: int = 180,
    fragment_sd: int = 20,
    error_rate: float = 0.002,
    stutter_rate: float = 0.0,
    seed: int = 13,
) -> dict[str, int]:
    registry = AlleleRegistry.load(registry_path)
    locus = registry.loci[pointer_id]
    weights: dict[int, float] = {}
    for part in weights_text.split(","):
        repeat_count, weight = part.split("=", 1)
        weights[int(repeat_count)] = float(weight)
    total = sum(weights.values())
    weights = {k: v / total for k, v in weights.items()}
    rng = random.Random(seed)
    counts = list(weights)
    probabilities = [weights[k] for k in counts]
    n_fragments = max(1, int(coverage * max(len(a.sequence) for a in locus.alleles) / max(1, 2 * read_length if out2 else read_length)))

    h1 = open_text(out1, "wt")
    h2 = open_text(out2, "wt") if out2 else None
    generated: dict[str, int] = {str(k): 0 for k in counts}
    try:
        for idx in range(n_fragments):
            repeat_count = rng.choices(counts, weights=probabilities, k=1)[0]
            observed_count = repeat_count
            if stutter_rate > 0 and rng.random() < stutter_rate:
                observed_count = max(1, repeat_count + rng.choice([-1, 1]))
            template = locus.left_flank + locus.motif * observed_count + locus.right_flank
            frag_len = min(len(template), max(read_length, int(rng.gauss(fragment_mean, fragment_sd))))
            start = rng.randint(0, max(0, len(template) - frag_len))
            fragment = template[start : start + frag_len]
            r1 = _mutate(fragment[:read_length], error_rate, rng)
            name = f"sim_{pointer_id}_{idx:06d}_L{repeat_count}"
            h1.write(f"@{name}/1\n{r1}\n+\n{'I'*len(r1)}\n")
            if h2:
                r2 = _mutate(reverse_complement(fragment[-read_length:]), error_rate, rng)
                h2.write(f"@{name}/2\n{r2}\n+\n{'I'*len(r2)}\n")
            generated[str(repeat_count)] += 1
    finally:
        h1.close()
        if h2:
            h2.close()
    return generated
