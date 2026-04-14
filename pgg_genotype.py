#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pgg_genotype_hybrid_fast.py
---------------------------------------------------------
Hybrid fast STR genotyper:
1. Keep the original fast banded SW stage for coarse scoring.
2. Run pHMM only on a small SW-shortlisted allele set.
3. Optionally cap very large candidate sets using frequency priors.
4. Early-cap per-locus FASTQ buffering to reduce memory.

This is meant to be much faster than the full pHMM version while still
using pHMM for final per-read allele likelihoods on the most plausible
candidates.
"""

import sys
import os
import gzip
import json
import math
import argparse
import csv
import random
import time
from collections import defaultdict

random.seed(1)
TRANS_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
NEG_INF = float("-inf")
LOG10_HALF = math.log10(0.5)


def revcomp(seq):
    return seq.translate(TRANS_TABLE)[::-1]


# -------------------------------------------------------
# Fast log helpers
# -------------------------------------------------------
def logaddexp_nat(a, b):
    if a == NEG_INF:
        return b
    if b == NEG_INF:
        return a
    if a < b:
        a, b = b, a
    return a + math.log1p(math.exp(b - a))


def logsumexp3_nat(a, b, c):
    return logaddexp_nat(logaddexp_nat(a, b), c)


def logsumexp4_nat(a, b, c, d):
    return logaddexp_nat(logaddexp_nat(a, b), logaddexp_nat(c, d))


def logsumexp_log10(values):
    vals = [v for v in values if v != NEG_INF and math.isfinite(v)]
    if not vals:
        return NEG_INF
    m = max(vals)
    return m + math.log10(sum(10 ** (v - m) for v in vals))


# -------------------------------------------------------
# Read / coordinate helpers
# -------------------------------------------------------
def normalize_read_id(raw_id):
    if raw_id.startswith("@"):
        raw_id = raw_id[1:]
    core = raw_id.split()[0]
    if core.endswith("/1") or core.endswith("/2"):
        core = core[:-2]
    elif core.endswith(".1") or core.endswith(".2"):
        core = core[:-2]
    return core


def parse_region(region_str):
    if not region_str:
        return None
    try:
        chrom, ran = region_str.split(":")
        start, end = map(int, ran.split("-"))
        return (chrom, start, end)
    except Exception:
        return None


def is_in_region(locus_id, region_tuple):
    if not region_tuple:
        return True
    try:
        parts = locus_id.split(":")
        l_chrom, l_start, l_end = None, 0, 0
        for p in parts:
            if "chr" in p:
                l_chrom = p
            elif "-" in p and p.replace("-", "").isdigit():
                s, e = map(int, p.split("-"))
                l_start, l_end = s, e
        if l_chrom is None or l_start == 0:
            return False
        r_chrom, r_start, r_end = region_tuple
        if l_chrom != r_chrom:
            return False
        return not (l_end < r_start or l_start > r_end)
    except Exception:
        return False


def parse_coords_from_gfa_id(locus_id):
    try:
        parts = locus_id.split(':')
        for i, p in enumerate(parts):
            if '-' in p and p.replace('-', '').isdigit():
                return f"{parts[i-1]}:{p}"
    except Exception:
        pass
    return None


def parse_coords_from_gaf_id(path_name):
    try:
        if path_name.startswith("allele_"):
            parts = path_name.split('_')
            if len(parts) >= 4:
                return f"{parts[1]}:{parts[2]}-{parts[3]}"
    except Exception:
        pass
    return None


# -------------------------------------------------------
# Prior helpers
# -------------------------------------------------------
def _cap_gq_from_delta(delta):
    if not math.isfinite(delta):
        return 99
    if delta < 0:
        delta = 0.0
    delta = min(delta, 9.9)
    return int(round(10.0 * delta))


def parse_popmix(s):
    if not s:
        return {}
    out = {}
    for kv in s.split(","):
        if not kv.strip():
            continue
        try:
            k, v = kv.split("=")
            out[k.strip().upper()] = float(v)
        except Exception:
            pass
    tot = sum(out.values())
    if tot > 0:
        for k in list(out.keys()):
            out[k] /= tot
    return out


def parse_extreme_list(lst):
    out = {}
    for s in lst or []:
        if "=" in s:
            lid, Ls = s.split("=", 1)
            try:
                out[lid.strip()] = int(Ls)
            except Exception:
                pass
    return out


def mix_prior_for_locus(freq_tab, locus_id, cand_L, popmix):
    eps = 1e-12
    raw = {}
    has_any = False
    locus_freqs = freq_tab.get(locus_id, {})
    for L in cand_L:
        pops = locus_freqs.get(L, {})
        if popmix:
            s = 0.0
            for p, w in popmix.items():
                s += w * pops.get(p, 0.0)
            raw[L] = s
        else:
            raw[L] = pops.get("ALL", 0.0)
        if raw[L] > 0:
            has_any = True

    if not has_any:
        for L in cand_L:
            raw[L] = 1.0

    tot = sum(raw.values()) + len(cand_L) * eps
    return {L: (v + eps) / tot for L, v in raw.items()}


def apply_extreme(prior, limit_val):
    eps = 1e-6
    n = len(prior)
    out = {}
    for L in prior:
        out[L] = 1.0 - (n - 1) * eps if L == limit_val else eps
    return out


def cap_candidates_by_prior(cands, prior, cap):
    if cap <= 0 or len(cands) <= cap:
        return sorted(cands)
    ranked = sorted(cands, key=lambda x: prior.get(x, 0.0), reverse=True)
    keep = ranked[:max(1, cap - 2)]
    keep.extend([min(cands), max(cands)])
    return sorted(set(keep))


# -------------------------------------------------------
# I/O
# -------------------------------------------------------
def parse_gaf_and_build_index(gaf_path, loci_map):
    print(f"[INFO] Loading GAF index from: {gaf_path}", file=sys.stderr)
    coord_to_gfa_id = {}
    for lid in loci_map:
        coord_key = parse_coords_from_gfa_id(lid)
        if coord_key:
            coord_to_gfa_id[coord_key] = lid

    mapping = defaultdict(set)
    op = gzip.open if str(gaf_path).endswith(".gz") else open
    try:
        with op(gaf_path, "rt", encoding="utf-8", errors="ignore") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 6:
                    continue
                raw_name = parts[0]
                path_name = parts[5]
                if path_name == "*" or not path_name:
                    continue

                target_lid = None
                if path_name in loci_map:
                    target_lid = path_name
                else:
                    c_key = parse_coords_from_gaf_id(path_name)
                    if c_key and c_key in coord_to_gfa_id:
                        target_lid = coord_to_gfa_id[c_key]

                if target_lid:
                    mapping[normalize_read_id(raw_name)].add(target_lid)
    except Exception as e:
        sys.exit(f"[ERROR] GAF Parsing Failed: {e}")
    print(f"[INFO] GAF Index ready. Found reads for {len(mapping)} unique IDs.", file=sys.stderr)
    return mapping


def stream_fastq_and_group(fq1, fq2, read_to_loci, per_locus_cap=None):
    loci_data = defaultdict(list)
    loci_counts = defaultdict(int)
    op1 = gzip.open if str(fq1).endswith(".gz") else open
    op2 = gzip.open if str(fq2).endswith(".gz") else open
    print("[INFO] Streaming FASTQ pairs...", file=sys.stderr)

    total_pairs = 0
    kept_pairs = 0

    try:
        with op1(fq1, "rt") as f1, op2(fq2, "rt") as f2:
            while True:
                h1 = f1.readline().strip()
                if not h1:
                    break
                s1 = f1.readline().strip()
                _ = f1.readline()
                q1 = f1.readline().strip()

                h2 = f2.readline().strip()
                s2 = f2.readline().strip()
                _ = f2.readline()
                q2 = f2.readline().strip()
                if not h2:
                    break

                total_pairs += 1
                if total_pairs % 5000000 == 0:
                    print(f"  Processed {total_pairs} pairs... (Kept {kept_pairs})", file=sys.stderr)

                core_id = normalize_read_id(h1)
                if core_id not in read_to_loci:
                    continue

                pair_data = ((h1, s1, q1), (h2, s2, q2))
                wrote_any = False
                for lid in read_to_loci[core_id]:
                    if per_locus_cap is not None and loci_counts[lid] >= per_locus_cap:
                        continue
                    loci_data[lid].append(pair_data)
                    loci_counts[lid] += 1
                    wrote_any = True
                if wrote_any:
                    kept_pairs += 1
    except Exception as e:
        sys.exit(f"[ERROR] FASTQ Reading Failed: {e}")

    print(f"[INFO] Extraction done. Kept {kept_pairs} relevant pairs.", file=sys.stderr)
    return loci_data


# -------------------------------------------------------
# Fast SW coarse stage
# -------------------------------------------------------
def sw_score(a, b, match=2, mismatch=-5, gap=-5, band=None):
    n, m = len(a), len(b)
    prev = [0] * (m + 1)
    best = 0

    for i in range(1, n + 1):
        cur = [0] * (m + 1)
        ai = a[i - 1]
        j_lo, j_hi = 1, m
        if band:
            j_lo = max(1, i - band)
            j_hi = min(m, i + band)
        for j in range(j_lo, j_hi + 1):
            sc = max(
                0,
                prev[j - 1] + (match if ai == b[j - 1] else mismatch),
                prev[j] + gap,
                cur[j - 1] + gap,
            )
            cur[j] = sc
            if sc > best:
                best = sc
        prev = cur
    return best


def coarse_sw_scores(seq, H_dict, args):
    if not seq:
        return None, None
    qseq = seq
    qseq_rc = revcomp(seq)
    scores = {}
    best = -1.0
    for L, hseq in H_dict.items():
        s = max(sw_score(qseq, hseq, band=args.band), sw_score(qseq_rc, hseq, band=args.band))
        scores[L] = s
        if s > best:
            best = s
    ideal_score = len(qseq) * 2
    if best < (ideal_score * args.min_score_ratio):
        return None, None
    return scores, best


def shortlist_from_sw(sw1, sw2, prior, args):
    order = []
    if sw1:
        order.extend(sorted(sw1, key=sw1.get, reverse=True)[:args.hybrid_top_k])
    if sw2:
        order.extend(sorted(sw2, key=sw2.get, reverse=True)[:args.hybrid_top_k])
    if prior and args.keep_top_prior > 0:
        order.extend(sorted(prior, key=prior.get, reverse=True)[:args.keep_top_prior])
    # Preserve order while deduplicating.
    seen = set()
    out = []
    for x in order:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def sw_log10_distribution(scores, args, allowed=None):
    if not scores:
        return None
    keys = allowed if allowed is not None else list(scores.keys())
    best = max(scores[L] for L in keys)
    tau = max(args.tau, 1e-6)
    vals = {L: math.exp((scores[L] - best) / tau) for L in keys}
    denom = sum(vals.values())
    out = {}
    for L in keys:
        p = max(vals[L] / denom, 1e-300)
        out[L] = math.log10(p)
    return out


# -------------------------------------------------------
# Optimized pHMM on shortlist only
# -------------------------------------------------------
def _clamp_prob(x):
    return min(max(x, 1e-12), 1.0 - 1e-12)


def build_allele_model(left, motif, repeat_count, right, args):
    motif = (motif or "N").upper()
    left = (left or "").upper()
    right = (right or "").upper()
    seq = left + motif * repeat_count + right
    m = len(seq)
    left_len = len(left)
    repeat_len = len(motif) * repeat_count
    repeat_start = left_len + 1
    repeat_end = left_len + repeat_len
    motif_len = max(1, len(motif))

    ref = " " + seq
    same_log = [0.0] * (m + 1)
    diff_log = [0.0] * (m + 1)
    m2m = [0.0] * (m + 1)
    m2i = [0.0] * (m + 1)
    m2d = [0.0] * (m + 1)
    i2m = [0.0] * (m + 1)
    i2i = [0.0] * (m + 1)
    d2m = [0.0] * (m + 1)
    d2d = [0.0] * (m + 1)

    flank_ins_emit = math.log(0.25)
    repeat_ins_emit = math.log(0.25)

    for j in range(1, m + 1):
        in_repeat = repeat_start <= j <= repeat_end
        boundary = in_repeat and ((j - repeat_start) % motif_len == 0)
        mismatch = args.repeat_mismatch if in_repeat else args.flank_mismatch
        same_log[j] = math.log(_clamp_prob(1.0 - mismatch))
        diff_log[j] = math.log(_clamp_prob(mismatch / 3.0))

        if in_repeat:
            m2i_p = args.repeat_m2i
            m2d_p = args.repeat_m2d
            i2i_p = args.repeat_i2i
            d2d_p = args.repeat_d2d
            if boundary:
                i2i_p = min(i2i_p + args.motif_boundary_bonus, 0.95)
                d2d_p = min(d2d_p + args.motif_boundary_bonus, 0.95)
        else:
            m2i_p = args.flank_m2i
            m2d_p = args.flank_m2d
            i2i_p = args.flank_i2i
            d2d_p = args.flank_d2d

        m2i_p = _clamp_prob(m2i_p)
        m2d_p = _clamp_prob(m2d_p)
        m2m_p = _clamp_prob(1.0 - m2i_p - m2d_p)
        i2i_p = _clamp_prob(i2i_p)
        i2m_p = _clamp_prob(1.0 - i2i_p)
        d2d_p = _clamp_prob(d2d_p)
        d2m_p = _clamp_prob(1.0 - d2d_p)

        m2m[j] = math.log(m2m_p)
        m2i[j] = math.log(m2i_p)
        m2d[j] = math.log(m2d_p)
        i2m[j] = math.log(i2m_p)
        i2i[j] = math.log(i2i_p)
        d2m[j] = math.log(d2m_p)
        d2d[j] = math.log(d2d_p)

    return {
        "seq": seq,
        "ref": ref,
        "m": m,
        "same_log": same_log,
        "diff_log": diff_log,
        "m2m": m2m,
        "m2i": m2i,
        "m2d": m2d,
        "i2m": i2m,
        "i2i": i2i,
        "d2m": d2m,
        "d2d": d2d,
        "flank_ins_emit": flank_ins_emit,
        "repeat_ins_emit": repeat_ins_emit,
        "repeat_start": repeat_start,
        "repeat_end": repeat_end,
    }


def phmm_forward_log10(seq, model):
    if not seq:
        return NEG_INF

    read = seq.upper()
    n = len(read)
    m = model["m"]
    ref = model["ref"]
    same_log = model["same_log"]
    diff_log = model["diff_log"]
    m2m = model["m2m"]
    m2i = model["m2i"]
    m2d = model["m2d"]
    i2m = model["i2m"]
    i2i = model["i2i"]
    d2m = model["d2m"]
    d2d = model["d2d"]
    rs = model["repeat_start"]
    re = model["repeat_end"]
    flank_ins_emit = model["flank_ins_emit"]
    repeat_ins_emit = model["repeat_ins_emit"]

    M_prev = [NEG_INF] * (m + 1)
    I_prev = [NEG_INF] * (m + 1)
    D_prev = [NEG_INF] * (m + 1)
    D_prev[0] = 0.0

    for j in range(1, m + 1):
        D_prev[j] = logsumexp3_nat(M_prev[j - 1] + m2d[j], D_prev[j - 1] + d2d[j], 0.0)

    for i in range(1, n + 1):
        base = read[i - 1]
        M_cur = [NEG_INF] * (m + 1)
        I_cur = [NEG_INF] * (m + 1)
        D_cur = [NEG_INF] * (m + 1)

        for j in range(1, m + 1):
            emit_M = same_log[j] if base == ref[j] else diff_log[j]
            emit_I = repeat_ins_emit if (rs <= j <= re) else flank_ins_emit

            M_cur[j] = emit_M + logsumexp4_nat(
                M_prev[j - 1] + m2m[j],
                I_prev[j - 1] + i2m[j],
                D_prev[j - 1] + d2m[j],
                0.0,
            )
            I_cur[j] = emit_I + logaddexp_nat(M_prev[j] + m2i[j], I_prev[j] + i2i[j])
            D_cur[j] = logaddexp_nat(M_cur[j - 1] + m2d[j], D_cur[j - 1] + d2d[j])

        M_prev, I_prev, D_prev = M_cur, I_cur, D_cur

    final_ln = NEG_INF
    for v in M_prev[1:]:
        final_ln = logaddexp_nat(final_ln, v)
    for v in I_prev[1:]:
        final_ln = logaddexp_nat(final_ln, v)
    for v in D_prev[1:]:
        final_ln = logaddexp_nat(final_ln, v)

    if final_ln == NEG_INF:
        return NEG_INF
    return final_ln / math.log(10.0)


def allele_read_log10_likelihood(seq, model, cache):
    key = (seq, id(model))
    if key in cache:
        return cache[key]
    fwd = phmm_forward_log10(seq, model)
    rev = phmm_forward_log10(revcomp(seq), model)
    out = logsumexp_log10([fwd + LOG10_HALF, rev + LOG10_HALF])
    cache[key] = out
    return out


def hybrid_read_log10_likelihood(seq, sw_scores, shortlist, allele_models, args, phmm_cache):
    if not seq or not sw_scores:
        return None

    # Stage 1: SW baseline over all candidates.
    sw_all_log10 = sw_log10_distribution(sw_scores, args)
    if sw_all_log10 is None:
        return None

    # Stage 2: pHMM only on shortlist.
    phmm_ll = {}
    best_ll = NEG_INF
    for L in shortlist:
        ll = allele_read_log10_likelihood(seq, allele_models[L], phmm_cache)
        phmm_ll[L] = ll
        if ll > best_ll:
            best_ll = ll

    if best_ll == NEG_INF:
        return sw_all_log10

    best_per_base = best_ll / max(len(seq), 1)
    if best_per_base < args.min_like_per_base:
        return sw_all_log10

    phmm_shifted = {L: (ll - best_ll) for L, ll in phmm_ll.items()}
    phmm_denom = logsumexp_log10(list(phmm_shifted.values()))
    phmm_short_probs = {L: 10 ** (v - phmm_denom) for L, v in phmm_shifted.items()}

    off = [L for L in sw_scores.keys() if L not in phmm_short_probs]
    off_mass = 0.0 if not off else min(max(args.off_shortlist_mass, 1e-8), 0.2)
    keep_mass = 1.0 - off_mass

    out = {}
    for L, p in phmm_short_probs.items():
        out[L] = math.log10(max(p * keep_mass, 1e-300))

    if off:
        sw_off_log10 = sw_log10_distribution(sw_scores, args, allowed=off)
        for L in off:
            p = off_mass * (10 ** sw_off_log10[L])
            out[L] = math.log10(max(p, 1e-300))

    return out


# -------------------------------------------------------
# GFA / frequency parsing
# -------------------------------------------------------
def parse_gfa_strs(path):
    loci = []
    try:
        with open(path, "r") as f:
            for ln in f:
                if not ln or ln.startswith("#"):
                    continue
                parts = ln.rstrip("\n").split("\t")
                if len(parts) >= 7 and parts[0] == "X" and parts[1] == "STR":
                    try:
                        cand_L = [int(x) for x in parts[6].split(",") if x]
                    except Exception:
                        cand_L = []
                    loci.append({
                        "locus_id": parts[2],
                        "motif": parts[3].upper(),
                        "left": parts[4].upper(),
                        "right": parts[5].upper(),
                        "candidates": cand_L,
                    })
    except Exception as e:
        sys.exit(f"[ERROR] GFA Parsing Failed: {e}")
    return loci


def load_freq_auto(path):
    tab = defaultdict(dict)
    if not path:
        return tab, "none"
    if not os.path.exists(path):
        print(f"[WARN] Frequency file not found: {path}, using uniform prior.", file=sys.stderr)
        return tab, "missing"
    try:
        with open(path, 'r') as f:
            for line in f:
                if not line.strip():
                    continue
                try:
                    obj = json.loads(line)
                    items = obj if isinstance(obj, list) else [obj]
                    for item in items:
                        lid = item.get("locus_id") or item.get("id") or item.get("pointer_id")
                        if not lid:
                            continue
                        data = item.get("L_freq") or item.get("freq") or item
                        for k, v in data.items():
                            if k.isdigit() and isinstance(v, dict):
                                tab[lid][int(k)] = v
                except Exception:
                    pass
        return tab, "json"
    except Exception:
        return tab, "error"


# -------------------------------------------------------
# Main
# -------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--gaf", required=True)
    ap.add_argument("--fq1", required=True)
    ap.add_argument("--fq2", required=True)
    ap.add_argument("--freq", required=False)
    ap.add_argument("--out", required=True)
    ap.add_argument("--popmix", required=False, default=None)
    ap.add_argument("--extreme", required=False, default=None, nargs="*")

    # Coarse SW stage (same interface as original script).
    ap.add_argument("--tau", type=float, default=12.0)
    ap.add_argument("--band", type=int, default=40)
    ap.add_argument("--target_good_reads", type=int, default=24)
    ap.add_argument("--max_reads_per_locus", type=int, default=48)
    ap.add_argument("--min_score_ratio", type=float, default=0.35)
    ap.add_argument("--region", required=False, help="Filter loci by region, e.g. chr19:40000000-50000000")

    # Hybrid speed controls.
    ap.add_argument("--hybrid_top_k", type=int, default=3, help="Run pHMM only on the top-k SW alleles per mate.")
    ap.add_argument("--keep_top_prior", type=int, default=1, help="Always include this many top-prior alleles in the pHMM shortlist.")
    ap.add_argument("--off_shortlist_mass", type=float, default=1e-3, help="Small probability mass reserved for non-shortlisted alleles.")
    ap.add_argument("--candidate_cap", type=int, default=16, help="Cap large locus candidate sets by prior; 0 disables.")
    ap.add_argument("--prefetch_pairs_per_locus", type=int, default=96, help="Early FASTQ buffering cap per locus.")

    # pHMM parameters.
    ap.add_argument("--flank_mismatch", type=float, default=0.02)
    ap.add_argument("--repeat_mismatch", type=float, default=0.05)
    ap.add_argument("--flank_m2i", type=float, default=0.03)
    ap.add_argument("--flank_m2d", type=float, default=0.03)
    ap.add_argument("--repeat_m2i", type=float, default=0.07)
    ap.add_argument("--repeat_m2d", type=float, default=0.07)
    ap.add_argument("--flank_i2i", type=float, default=0.30)
    ap.add_argument("--flank_d2d", type=float, default=0.30)
    ap.add_argument("--repeat_i2i", type=float, default=0.55)
    ap.add_argument("--repeat_d2d", type=float, default=0.55)
    ap.add_argument("--motif_boundary_bonus", type=float, default=0.08)
    ap.add_argument("--min_like_per_base", type=float, default=-0.90)

    args = ap.parse_args()

    out_dir = os.path.dirname(args.out)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    popmix = parse_popmix(args.popmix)
    extreme_map = parse_extreme_list(args.extreme)
    region_filter = parse_region(args.region)
    if region_filter:
        print(f"[INFO] Region Filter Enabled: {args.region}", file=sys.stderr)

    all_loci = parse_gfa_strs(args.gfa)
    loci_map = {}
    for l in all_loci:
        if is_in_region(l["locus_id"], region_filter):
            loci_map[l["locus_id"]] = l
    print(f"[INFO] Loci active after filtering: {len(loci_map)} (Original GFA: {len(all_loci)})", file=sys.stderr)
    if not loci_map:
        sys.exit("[ERROR] No loci found in the specified region! Check your --region or GFA.")

    freq_tab, _ftype = load_freq_auto(args.freq)
    read_whitelist = parse_gaf_and_build_index(args.gaf, loci_map)
    grouped_reads = stream_fastq_and_group(
        args.fq1,
        args.fq2,
        read_whitelist,
        per_locus_cap=max(args.prefetch_pairs_per_locus, args.max_reads_per_locus),
    )

    out_rows = []
    total_active_loci = len(grouped_reads)
    print(f"[INFO] Genotyping {total_active_loci} active loci with reads...", file=sys.stderr)

    start_time = time.time()
    processed_count = 0

    for lid in loci_map:
        if lid not in grouped_reads:
            continue

        processed_count += 1
        if processed_count % 50 == 0:
            elapsed = time.time() - start_time
            print(f"  [Progress] {processed_count}/{total_active_loci} loci processed ({elapsed:.1f}s)...", file=sys.stderr)

        pairs_for_this_locus = grouped_reads[lid]
        loc = loci_map[lid]
        cands = sorted(set(loc["candidates"]))
        if not cands:
            cands = list(range(3, 20))

        prior = mix_prior_for_locus(freq_tab, lid, cands, popmix)
        if lid in extreme_map:
            prior = apply_extreme(prior, extreme_map[lid])
        cands = cap_candidates_by_prior(cands, prior, args.candidate_cap)
        prior = {L: prior[L] for L in cands}
        # renormalize after candidate capping
        tot_prior = sum(prior.values())
        if tot_prior > 0:
            prior = {L: v / tot_prior for L, v in prior.items()}
        else:
            prior = {L: 1.0 / len(cands) for L in cands}

        H = {L: (loc["left"] + loc["motif"] * L + loc["right"]) for L in cands}
        allele_models = {L: build_allele_model(loc["left"], loc["motif"], L, loc["right"], args) for L in cands}
        phmm_cache = {}

        per_pair_joint_log10P = []
        valid_pairs_count = 0
        total_checked = 0

        for (r1_data, r2_data) in pairs_for_this_locus:
            if valid_pairs_count >= args.target_good_reads:
                break
            if total_checked >= args.max_reads_per_locus:
                break
            total_checked += 1

            (_n1, s1, _q1) = r1_data
            (_n2, s2, _q2) = r2_data

            sw1, _best1 = coarse_sw_scores(s1, H, args)
            sw2, _best2 = coarse_sw_scores(s2, H, args)
            if sw1 is None and sw2 is None:
                continue

            shortlist = shortlist_from_sw(sw1, sw2, prior, args)
            if not shortlist:
                continue

            logP_R1 = hybrid_read_log10_likelihood(s1, sw1, shortlist, allele_models, args, phmm_cache) if sw1 else None
            logP_R2 = hybrid_read_log10_likelihood(s2, sw2, shortlist, allele_models, args, phmm_cache) if sw2 else None
            if logP_R1 is None and logP_R2 is None:
                continue

            valid_pairs_count += 1
            joint_logP = {}
            for L in cands:
                s = 0.0
                if logP_R1 is not None:
                    s += logP_R1[L]
                if logP_R2 is not None:
                    s += logP_R2[L]
                joint_logP[L] = s
            per_pair_joint_log10P.append(joint_logP)

        if valid_pairs_count == 0:
            out_rows.append({"locus_id": lid, "GT": "./.", "GQ": 0, "n_reads": 0})
            continue

        post = {}
        cand_pairs = []
        for i, L1 in enumerate(cands):
            for L2 in cands[i:]:
                cand_pairs.append((L1, L2))

        for (L1, L2) in cand_pairs:
            ll = 0.0
            for d in per_pair_joint_log10P:
                p_l1 = 10 ** d[L1]
                p_l2 = 10 ** d[L2]
                pr = 0.5 * p_l1 + 0.5 * p_l2
                ll += math.log10(max(pr, 1e-300))
            pg = (prior[L1] ** 2) if L1 == L2 else (2 * prior[L1] * prior[L2])
            post[(L1, L2)] = ll + math.log10(max(pg, 1e-300))

        sorted_post = sorted(post.items(), key=lambda x: x[1], reverse=True)
        (gtL1, gtL2), best_p = sorted_post[0]
        second_p = sorted_post[1][1] if len(sorted_post) > 1 else -999.0
        GQ = _cap_gq_from_delta(best_p - second_p)

        out_rows.append({
            "locus_id": lid,
            "GT": f"{gtL1}/{gtL2}",
            "GQ": GQ,
            "n_reads": valid_pairs_count,
        })

    keys = ["locus_id", "GT", "GQ", "n_reads"]
    try:
        with open(args.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
            w.writeheader()
            for r in out_rows:
                w.writerow(r)
        print(f"[INFO] Success. Output written to {args.out}", file=sys.stderr)
    except Exception as e:
        sys.exit(f"[ERROR] Writing output failed: {e}")


if __name__ == "__main__":
    main()
