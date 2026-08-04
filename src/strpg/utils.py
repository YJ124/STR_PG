from __future__ import annotations

import gzip
import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Iterator, Sequence, TextIO, TypeVar

import numpy as np

T = TypeVar("T")

DNA_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]


def open_text(path: str | Path, mode: str = "rt") -> TextIO:
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, mode)  # type: ignore[return-value]
    return open(p, mode, encoding=None if "b" in mode else "utf-8")  # type: ignore[return-value]


@lru_cache(maxsize=1 << 16)
def stable_hash(text: str, digest_size: int = 8) -> int:
    return int.from_bytes(hashlib.blake2b(text.encode("ascii"), digest_size=digest_size).digest(), "big")


def canonical_kmer(kmer: str) -> tuple[str, int]:
    kmer = kmer.upper()
    rc = reverse_complement(kmer)
    if rc < kmer:
        return rc, -1
    return kmer, 1


def logsumexp(values: Sequence[float] | np.ndarray) -> float:
    if len(values) == 0:
        return -math.inf
    arr = np.asarray(values, dtype=float)
    m = float(np.max(arr))
    if not math.isfinite(m):
        return m
    return m + math.log(float(np.exp(arr - m).sum()))


def log_normalize(log_values: Sequence[float]) -> list[float]:
    z = logsumexp(log_values)
    if not math.isfinite(z):
        n = len(log_values)
        return [1.0 / n] * n if n else []
    return [math.exp(v - z) for v in log_values]


def atomic_json_dump(obj: object, path: str | Path, indent: int = 2) -> None:
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(obj, handle, ensure_ascii=False, indent=indent, sort_keys=True)
        handle.write("\n")
    tmp.replace(path)


def chunks(items: Sequence[T], size: int) -> Iterator[Sequence[T]]:
    for i in range(0, len(items), size):
        yield items[i : i + size]


def phred_from_error_probability(error_probability: float, cap: int = 99) -> int:
    p = min(max(error_probability, 1e-300), 1.0)
    return min(cap, int(round(-10.0 * math.log10(p))))


def safe_float(value: str | float | int | None, default: float = 0.0) -> float:
    try:
        return float(value) if value is not None else default
    except (TypeError, ValueError):
        return default


def ensure_parent(path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
