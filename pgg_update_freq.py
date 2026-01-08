#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pgg_update_freq.py

Incrementally updates the population frequency database `freq.jsonl`
based on the STR genotyping results (genotypes.tsv) of a single sample
(using homozygous loci only).

Logic:
  - Only count homozygous genotypes (GT is x/x).
  - Each homozygous sample contributes 2 copies of that allele.
  - Update `counts[pop][allele] += 2` and `total_alleles[pop] += 2`
    for the corresponding `pointer_id` + population label in `freq.jsonl`.
    Then recalculate:
        freq[pop][allele] = counts[pop][allele] / total_alleles[pop]

Notes:
  - To be compatible with existing `freq.jsonl`, the script will:
      * If the original record lacks `counts` or `total_alleles` fields,
        generate "pseudo-counts" using `freq * init_total_alleles`
        before proceeding with the incremental update.
  - `locus_id` and `pointer_id` are assumed to be the same string by default
    (consistent with your existing pipeline).

Usage Example:

  python pgg_update_freq.py \\
      --geno genotypes.tsv \\
      --freq-in freq19.jsonl \\
      --freq-out freq19.updated.jsonl \\
      --pop EAS \\
      --min-GQ 20 \\
      --init-total-alleles 1000
"""

import os
import json
import argparse
from typing import Dict, Any, Tuple


def load_freq_jsonl(path: str,
                    init_total_alleles: int) -> Dict[str, Dict[str, Any]]:
    """
    Reads `freq.jsonl` and returns:
        { pointer_id: { 'pointer_id':..., 'freq':..., 'counts':..., 'total_alleles':... } }
    If `counts` or `total_alleles` are missing in the original jsonl,
    generate pseudo-counts using `freq * init_total_alleles`.
    """
    db: Dict[str, Dict[str, Any]] = {}
    if not os.path.exists(path):
        return db

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            pid = obj.get("pointer_id")
            if not pid:
                continue

            freq = obj.get("freq", {})
            counts = obj.get("counts", {})
            total_alleles = obj.get("total_alleles", {})

            # If counts/total_alleles are missing, derive an initial count from freq
            if not counts or not total_alleles:
                counts = {}
                total_alleles = {}
                for pop, dist in freq.items():
                    pop_counts = {}
                    for L_str, p in dist.items():
                        # Use frequency * init_total_alleles to get an approximate integer count
                        c = max(0, int(round(float(p) * init_total_alleles)))
                        if c > 0:
                            pop_counts[L_str] = c
                    total_c = sum(pop_counts.values())
                    # If a population has no allele counts, total_alleles will be 0
                    counts[pop] = pop_counts
                    total_alleles[pop] = total_c

            obj["freq"] = freq
            obj["counts"] = counts
            obj["total_alleles"] = total_alleles
            db[pid] = obj

    return db


def save_freq_jsonl(db: Dict[str, Dict[str, Any]], out_path: str) -> None:
    """
    Writes the updated frequency database back to a JSONL file.
    Preserves `pointer_id`, `freq`, `counts`, and `total_alleles` fields.
    """
    with open(out_path, "w", encoding="utf-8") as f:
        for pid, obj in db.items():
            # Ensure pointer_id consistency
            obj["pointer_id"] = pid
            # Recalculate freq based on counts & total_alleles first
            freq = obj.get("freq", {})
            counts = obj.get("counts", {})
            total_alleles = obj.get("total_alleles", {})

            new_freq = {}
            for pop, alle_counts in counts.items():
                ta = total_alleles.get(pop, 0)
                if ta <= 0:
                    continue
                pop_freq = {}
                for L_str, c in alle_counts.items():
                    if c <= 0:
                        continue
                    pop_freq[L_str] = c / float(ta)
                if pop_freq:
                    new_freq[pop] = pop_freq
            obj["freq"] = new_freq

            f.write(json.dumps(obj, ensure_ascii=False) + "\n")


def parse_gt(gt: str) -> Tuple[int, int]:
    """
    Parses a GT string like '10/12' into (10, 12).
    If it is './.' or invalid, raises an exception to be handled by the caller.
    """
    if gt is None:
        raise ValueError("GT is None")
    gt = gt.strip()
    if gt == "./.":
        raise ValueError("missing GT")
    if "/" not in gt:
        raise ValueError(f"bad GT format: {gt}")
    a, b = gt.split("/", 1)
    return int(a), int(b)


def update_freq_with_sample(
    db: Dict[str, Dict[str, Any]],
    geno_tsv: str,
    pop: str,
    min_gq: int = 0
) -> None:
    """
    Core function: Updates the frequency database `db` using a single sample's `genotypes.tsv`.
    Only counts:
      - Homozygous GT (L/L)
      - GQ >= min_gq
    """
    if pop is None or pop == "":
        raise ValueError("pop (population label) cannot be empty, e.g., EAS/EUR/AFR.")

    with open(geno_tsv, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        # Expected column names (consistent with your original pgg_genotype.tsv):
        # locus_id, motif, n_reads, candidates, GT, GQ, post_best_log10, post2_log10, APL_json
        col_idx = {name: i for i, name in enumerate(header)}

        required_cols = ["locus_id", "GT", "GQ"]
        for col in required_cols:
            if col not in col_idx:
                raise RuntimeError(
                    f"Missing required column in TSV: {col}, current header: {header}"
                )

        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            locus_id = parts[col_idx["locus_id"]]
            gt_str = parts[col_idx["GT"]]
            gq_str = parts[col_idx["GQ"]]

            # Skip missing GT
            if gt_str == "./.":
                continue

            # GQ filtering
            try:
                gq = int(gq_str)
            except Exception:
                gq = 0
            if gq < min_gq:
                continue

            # Parse GT
            try:
                L1, L2 = parse_gt(gt_str)
            except Exception:
                continue

            # Only count homozygous loci
            if L1 != L2:
                continue
            L = L1

            # Start incremental update on the frequency database
            pid = locus_id  # Default: locus_id == pointer_id
            if pid not in db:
                # If the locus doesn't exist in the original freq.jsonl, create a new record
                db[pid] = {
                    "pointer_id": pid,
                    "freq": {},
                    "counts": {},
                    "total_alleles": {}
                }

            obj = db[pid]
            freq = obj.setdefault("freq", {})
            counts = obj.setdefault("counts", {})
            total_alleles = obj.setdefault("total_alleles", {})

            pop_counts = counts.setdefault(pop, {})
            pop_total = total_alleles.get(pop, 0)

            L_str = str(L)

            # For homozygous individuals: contribute 2 copies of this allele
            c_old = pop_counts.get(L_str, 0)
            c_new = c_old + 2
            pop_counts[L_str] = c_new
            pop_total += 2
            total_alleles[pop] = pop_total

            # Frequency is not calculated here; only update counts / total_alleles.
            # Frequency will be recalculated when writing out the file.
            obj["counts"] = counts
            obj["total_alleles"] = total_alleles
            db[pid] = obj


def main():
    ap = argparse.ArgumentParser(
        description="Incrementally update freq.jsonl based on homozygous loci in genotypes.tsv"
    )
    ap.add_argument("--geno", required=True,
                    help="Path to genotypes.tsv output from pgg_genotype")
    ap.add_argument("--freq-in", required=True,
                    help="Path to original freq.jsonl / freq19.jsonl")
    ap.add_argument("--freq-out", required=True,
                    help="Path to output updated freq.jsonl")
    ap.add_argument("--pop", required=True,
                    help="Population label for this sample (e.g., EAS/EUR/AFR)")
    ap.add_argument("--min-GQ", type=int, default=0,
                    help="Only use loci with GQ >= this threshold (default 0)")
    ap.add_argument("--init-total-alleles", type=int, default=1000,
                    help="Base count for generating pseudo-counts when original freq.jsonl "
                         "lacks counts/total_alleles (default 1000)")

    args = ap.parse_args()

    # 1) Load old frequency database
    print(f"[INFO] loading freq jsonl from: {args.freq_in}")
    db = load_freq_jsonl(args.freq_in, init_total_alleles=args.init_total_alleles)
    print(f"[INFO] loaded {len(db)} loci from freq jsonl.")

    # 2) Update frequency with this sample's genotypes.tsv
    print(f"[INFO] updating freq with genotypes from: {args.geno}")
    update_freq_with_sample(
        db,
        geno_tsv=args.geno,
        pop=args.pop,
        min_gq=args.min_GQ
    )

    # 3) Write out new freq.jsonl
    print(f"[INFO] writing updated freq jsonl to: {args.freq_out}")
    save_freq_jsonl(db, args.freq_out)
    print("[DONE] freq jsonl updated.")


if __name__ == "__main__":
    main()
