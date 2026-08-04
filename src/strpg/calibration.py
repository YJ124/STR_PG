from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from .utils import open_text


def calibrate_threshold(
    genotype_tsv: str | Path,
    out_json: str | Path,
    column: str = "D_PAPER",
    quantile: float = 0.99,
) -> dict:
    if not (0.0 < quantile < 1.0):
        raise ValueError("quantile must be in (0,1)")
    values: list[float] = []
    with open_text(genotype_tsv, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if column not in (reader.fieldnames or []):
            raise ValueError(f"Column {column} not found")
        for row in reader:
            try:
                values.append(float(row[column]))
            except (TypeError, ValueError):
                continue
    if not values:
        raise ValueError("No numeric calibration scores found")
    threshold = float(np.quantile(np.asarray(values), quantile, method="higher"))
    payload = {
        "metric": column,
        "quantile": quantile,
        "threshold": threshold,
        "n_null_loci": len(values),
        "note": "Empirical null threshold; use a null dataset matched for coverage and candidate-set size.",
    }
    with open(out_json, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return payload


def load_threshold(path: str | Path) -> float:
    with open(path, "r", encoding="utf-8") as handle:
        return float(json.load(handle)["threshold"])
