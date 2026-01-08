#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pgg_genotype_str.py  (JSONL frequency support version, supports single/paired-end FASTQ)
---------------------------------------------------------
A simplified, runnable small-scale STR genotyper (diploid), supporting:
- webSTR style TSV frequency tables (locus_id, L, AFR/EUR/.../ALL)
- JSONL/JSON frequency tables (various common nested structures, see load_freq_auto comments)
- Admixture priors (--popmix "EUR=0.6,AFR=0.4")
- Extreme priors (--extreme "locus_id=L")
- Single-end FASTQ (--fq) or Paired-end FASTQ (--fq1 + --fq2)

Output: One line per locus in TSV format, containing GT/GQ, like/post scores, prior debug columns, APL_json.
This script is used to "verify if priors are entering the scoring pipeline", not a production-optimized version.
"""
import sys, os, gzip, json, math, argparse, csv, random
from collections import defaultdict

random.seed(1)

def _cap_gq_from_delta(delta):
    """Convert Δ=post_best-post2 to GQ; safe for inf/NaN, and capped at 99."""
    if not math.isfinite(delta):
        return 99
    if delta < 0:
        delta = 0.0
    # Set an upper limit to avoid excessive values (10 * 9.9 = 99)
    delta = min(delta, 9.9)
    return int(round(10.0 * delta))

# ------------------ FASTQ Reading ------------------
def read_fastq(path):
    """Read single FASTQ(.gz): yield (name, seq, qual)"""
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        while True:
            name = f.readline().rstrip()
            if not name:
                break
            seq  = f.readline().rstrip()
            plus = f.readline().rstrip()
            qual = f.readline().rstrip()
            if not qual:
                break
            yield name[1:] if name.startswith("@") else name, seq, qual

def load_reads_from_args(args):
    """
    Load reads based on command line arguments:
      - Single-end: --fq
      - Paired-end: --fq1 + --fq2

    In paired-end mode, every read in R1 and R2 is treated as an independent read for SW scoring.
    To distinguish them, /1 or /2 is appended to the name (if not already present).
    """
    if args.fq and (args.fq1 or args.fq2):
        print("[ERROR] Cannot use both --fq and --fq1/--fq2, please choose one.", file=sys.stderr)
        sys.exit(1)
    if (args.fq1 and not args.fq2) or (args.fq2 and not args.fq1):
        print("[ERROR] Paired-end mode requires both --fq1 and --fq2.", file=sys.stderr)
        sys.exit(1)
    if not args.fq and not (args.fq1 and args.fq2):
        print("[ERROR] Please use either --fq (single-end) or --fq1 + --fq2 (paired-end).", file=sys.stderr)
        sys.exit(1)

    reads = []
    if args.fq:
        # Single-end
        print(f"[INFO] Using single-end FASTQ: {args.fq}", file=sys.stderr)
        for name, seq, qual in read_fastq(args.fq):
            reads.append((name, seq, qual))
    else:
        # Paired-end
        print(f"[INFO] Using paired-end FASTQ: R1={args.fq1}  R2={args.fq2}", file=sys.stderr)
        # R1
        for name, seq, qual in read_fastq(args.fq1):
            if not name.endswith("/1"):
                name = name + "/1"
            reads.append((name, seq, qual))
        # R2
        for name, seq, qual in read_fastq(args.fq2):
            if not name.endswith("/2"):
                name = name + "/2"
            reads.append((name, seq, qual))

    if not reads:
        print("[ERROR] No reads loaded from FASTQ, is the file empty?", file=sys.stderr)
        sys.exit(2)
    print(f"[INFO] Total reads loaded: {len(reads)}", file=sys.stderr)
    return reads

# ------------------ GFA Parsing (X STR lines) ------------------
def parse_gfa_strs(path):
    loci = []  # list of dict
    with open(path, "r") as f:
        for ln in f:
            if not ln or ln[0] == "#":
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) >= 7 and parts[0] == "X" and parts[1] == "STR":
                locus_id, motif, left, right, LsCSV = parts[2], parts[3], parts[4], parts[5], parts[6]
                try:
                    cand_L = [int(x) for x in LsCSV.split(",") if x]
                except:
                    cand_L = []
                loci.append({
                    "locus_id": locus_id,
                    "motif": motif.upper(),
                    "left": left.upper(),
                    "right": right.upper(),
                    "candidates": cand_L,
                })
    return loci

# ------------------ Frequency Reading (Auto-detect TSV / JSONL / JSON) ------------------
POPS = ["AFR","EUR","EAS","AMR","SAS","ALL"]

import csv, json

def load_freq_tsv(path):
    tab = defaultdict(dict)  # tab[locus_id][L] = {pop:val}
    with open(path, newline="", encoding="utf-8") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            lid = (row.get("locus_id") or row.get("id") or row.get("locus") or "").strip()
            if not lid:
                continue
            try:
                L = int(row.get("L"))
            except:
                continue
            d = {}
            for p in POPS:
                v = row.get(p, "")
                if v not in (None, "", "NA"):
                    try:
                        d[p] = float(v)
                    except:
                        pass
            if d:
                tab[lid][L] = d
    return tab

def _ingest_row_like(tab, obj):
    """Extract frequencies from a single JSON object 'obj' as robustly as possible.
    Supports:
      A) Flat: {"locus_id": "...", "L": 10, "ALL": 0.1, "EUR": 0.2, ...}
      B) pop->L: {"locus_id":"...", "freq": {"ALL":{"10":0.1,"11":0.2}, "EUR":{...}}}
      C) L->pop: {"locus_id":"...", "L_freq": {"10":{"ALL":0.1,"EUR":0.2}, "11":{...}}}
      D) Container key aliases: "priors"/"allele_freqs"/"pops"
    """
    lid = (obj.get("locus_id") or obj.get("pointer_id") or obj.get("id") or obj.get("locus") or "").strip()
    if not lid:
        return

    # A) Flat
    if "L" in obj and any(k in obj for k in POPS):
        try:
            L = int(obj["L"])
            d = {}
            for p in POPS:
                if p in obj:
                    try:
                        d[p] = float(obj[p])
                    except:
                        pass
            if d:
                tab[lid][L] = d
                return
        except:
            pass

    # B/C/D) Container
    cand_container_keys = ["freq","L_freq","priors","allele_freqs","pops"]
    container = None
    for k in cand_container_keys:
        if k in obj and isinstance(obj[k], dict):
            container = obj[k]
            break
    if container is None:
        # Could also be that pop keys are directly {L:val}
        has_pop_key = any((k in obj and isinstance(obj[k], dict)) for k in POPS)
        if has_pop_key:
            for p in POPS:
                if p in obj and isinstance(obj[p], dict):
                    for Ls, val in obj[p].items():
                        try:
                            L = int(Ls); val = float(val)
                            tab[lid].setdefault(L, {})[p] = val
                        except:
                            pass
            return
        # Or obj itself is {L:{pop:val}}
        try:
            for Ls, mp in obj.items():
                L = int(Ls)
                if isinstance(mp, dict):
                    for p, val in mp.items():
                        if p in POPS:
                            tab[lid].setdefault(L, {})[p] = float(val)
            return
        except:
            return
    else:
        # container could be pop->L or L->pop
        if any(k in POPS for k in container.keys()):
            # pop->L
            for p, mp in container.items():
                if p not in POPS:
                    continue
                if isinstance(mp, dict):
                    it = mp.items()
                elif isinstance(mp, list):
                    it = []
                    for item in mp:
                        if isinstance(item, dict) and "L" in item:
                            it.append((item["L"], item.get("v", item.get("val", item.get("freq", None)))))
                else:
                    continue
                for Ls, val in it:
                    try:
                        L = int(Ls); val = float(val)
                        tab[lid].setdefault(L, {})[p] = val
                    except:
                        pass
        else:
            # L->pop
            for Ls, mp in container.items():
                try:
                    L = int(Ls)
                except:
                    continue
                if not isinstance(mp, dict):
                    continue
                for p, val in mp.items():
                    if p in POPS:
                        try:
                            tab[lid].setdefault(L, {})[p] = float(val)
                        except:
                            pass

def load_freq_json_lines(path, force_jsonl=False):
    tab = defaultdict(dict)

    if force_jsonl:
        # Strict JSONL: parse line by line
        with open(path, "r", encoding="utf-8") as f:
            for ln in f:
                ln = ln.strip()
                if not ln or ln.startswith("#"):
                    continue
                try:
                    obj = json.loads(ln)
                    if isinstance(obj, (dict, list)):
                        if isinstance(obj, list):
                            for x in obj:
                                if isinstance(x, dict):
                                    _ingest_row_like(tab, x)
                        else:
                            _ingest_row_like(tab, obj)
                except:
                    continue
        return tab

    # AUTO: try whole JSON first; fail then fallback to JSONL
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)   # If it's JSONL, this will likely raise Extra data
        # Success: data might be list or dict
        if isinstance(data, list):
            for x in data:
                if isinstance(x, dict):
                    _ingest_row_like(tab, x)
        elif isinstance(data, dict):
            _ingest_row_like(tab, data)
        return tab
    except Exception:
        # Fallback JSONL
        return load_freq_json_lines(path, force_jsonl=True)

def load_freq_auto(path):
    low = path.lower()
    if low.endswith(".jsonl"):
        return load_freq_json_lines(path, force_jsonl=True), "json"
    elif low.endswith(".json"):
        return load_freq_json_lines(path, force_jsonl=False), "json"
    else:
        return load_freq_tsv(path), "tsv"

# ------------------ Admixture and Extreme Priors ------------------
def parse_popmix(s):
    # "EUR=0.6,AFR=0.4" -> dict
    if not s:
        return {}
    out = {}
    for kv in s.split(","):
        if not kv.strip():
            continue
        k,v = kv.split("=")
        out[k.strip().upper()] = float(v)
    # Normalize
    tot = sum(out.values())
    if tot > 0:
        for k in out:
            out[k] /= tot
    return out

def mix_prior_for_locus(freq_tab, locus_id, cand_L, popmix):
    # Returns p_L (smoothed, normalized)
    eps = 1e-12
    raw = {}
    has_any = False
    for L in cand_L:
        pops = freq_tab.get(locus_id, {}).get(L, {})
        if popmix:
            s = 0.0
            for p,w in popmix.items():
                s += w * pops.get(p, 0.0)
            raw[L] = s
        else:
            # Default ALL
            raw[L] = pops.get("ALL", 0.0)
        has_any = has_any or (raw[L] > 0)
    # Smoothing: if all are 0, revert to uniform
    if not has_any:
        for L in cand_L:
            raw[L] = 1.0
    # Add eps to prevent zero
    for L in cand_L:
        raw[L] = raw[L] + eps
    # Normalize
    tot = sum(raw.values())
    for L in cand_L:
        raw[L] /= tot
    return raw

def parse_extreme_list(lst):
    # ["locus_id=L", "locus2=7", ...]
    out = {}
    for s in lst or []:
        if "=" in s:
            lid, Ls = s.split("=",1)
            try:
                out[lid] = int(Ls)
            except:
                pass
    return out

def apply_extreme(prior_dict, extreme_L):
    # Approximating prior change to: L*=1-ε, others=(ε/(n-1))
    if extreme_L is None:
        return prior_dict
    ks = list(prior_dict.keys())
    n = len(ks)
    if extreme_L not in ks or n<=1:
        return prior_dict
    eps = 1e-6
    for L in ks:
        prior_dict[L] = eps/(n-1)
    prior_dict[extreme_L] = 1.0 - eps
    return prior_dict

# ------------------ Simple Dual-Ended SW (linear gap; only calculate max score) ------------------
def sw_score(a, b, match=2, mismatch=-3, gap=-4, band=None):
    # Returns max local alignment score (no traceback)
    # Optional band: if given, limits |i - j| <= band (saves time)
    n, m = len(a), len(b)
    prev = [0]*(m+1)
    best = 0
    for i in range(1, n+1):
        cur = [0]*(m+1)
        ai = a[i-1]
        j_lo = 1
        j_hi = m
        if band is not None:
            j_lo = max(1, i - band)
            j_hi = min(m, i + band)
        for j in range(1, j_lo):
            cur[j] = 0
        for j in range(j_lo, j_hi+1):
            sc_diag = prev[j-1] + (match if ai == b[j-1] else mismatch)
            sc_up   = prev[j] + gap
            sc_left = cur[j-1] + gap
            sc = max(0, sc_diag, sc_up, sc_left)
            cur[j] = sc
            if sc > best:
                best = sc
        for j in range(j_hi+1, m+1):
            cur[j] = 0
        prev = cur
    return best

# ------------------ Genotype Traversal and Scoring ------------------
def all_genotype_pairs(cands):
    out = []
    for i, L1 in enumerate(cands):
        for j, L2 in enumerate(cands[i:]):
            out.append((L1, L2))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True, help="GFA file containing X STR lines")
    # Single-end or Paired-end, choose one:
    ap.add_argument("--fq",  required=False, help="Single-end FASTQ(.gz)")
    ap.add_argument("--fq1", required=False, help="Paired-end R1 FASTQ(.gz)")
    ap.add_argument("--fq2", required=False, help="Paired-end R2 FASTQ(.gz)")
    ap.add_argument("--freq", required=False, default=None, help="Allele frequency table, supports .tsv / .jsonl / .json")
    ap.add_argument("--popmix", required=False, default=None, help='e.g. \"EUR=0.6,AFR=0.4\", otherwise uses ALL column or uniform')
    ap.add_argument("--extreme", required=False, default=None, nargs="*", help='List in format locus_id=L, for extreme prior verification')
    ap.add_argument("--tau", type=float, default=12.0, help="Softmax temperature (scaling SW score to probability)")
    ap.add_argument("--band", type=int, default=60, help="SW bandwidth, limits |i-j| to speed up")
    ap.add_argument("--max_reads_per_locus", type=int, default=500, help="Max reads to use per locus (random subsampling)")
    ap.add_argument("--out", required=True, help="Output TSV file path")
    args = ap.parse_args()

    loci = parse_gfa_strs(args.gfa)
    if not loci:
        print("No X STR lines found in GFA", file=sys.stderr); sys.exit(2)

    freq_tab = {}
    freq_type = None
    if args.freq:
        freq_tab, freq_type = load_freq_auto(args.freq)
    popmix = parse_popmix(args.popmix)
    extreme_map = parse_extreme_list(args.extreme)

    # Load all reads into memory (small samples are fine; paired-end means R1+R2)
    reads = load_reads_from_args(args)

    out_rows = []
    for loc in loci:
        lid   = loc["locus_id"]
        motif = loc["motif"]
        left  = loc["left"]
        right = loc["right"]
        cands = sorted(set(loc["candidates"])) or list(range(3, 20))

        # H_L Sequence Pool
        H = {L: (left + motif*L + right) for L in cands}

        # Sample reads
        idxs = list(range(len(reads)))
        random.shuffle(idxs)
        idxs = idxs[:args.max_reads_per_locus]
        reads_sel = [reads[i] for i in idxs]

        # Per read: Softmax SW scores into P(r|L)
        per_read_log10PL = []
        for (qname, qseq, qqual) in reads_sel:
            # Calculate SW score for each L
            scores = {}
            bestS = None
            for L, hseq in H.items():
                s = sw_score(qseq, hseq, match=2, mismatch=-3, gap=-4, band=args.band)
                scores[L] = s
                if (bestS is None) or (s > bestS):
                    bestS = s
            # softmax -> probability (saved as log10)
            log10P = {}
            denom = 0.0
            for L, s in scores.items():
                val = math.exp((s - bestS)/max(args.tau, 1e-6))
                log10P[L] = val
                denom += val
            for L in scores.keys():
                p = max(log10P[L]/denom, 1e-300)
                log10P[L] = math.log10(p)
            per_read_log10PL.append(log10P)

        n_reads = len(per_read_log10PL)

        # APL_json (Aggregated "cost" per allele): -10 * (sum_r log10 P(r|L) - max_L ...)
        sum_log10_P_L = {L: 0.0 for L in cands}
        for d in per_read_log10PL:
            for L,v in d.items():
                sum_log10_P_L[L] += v
        max_sum = max(sum_log10_P_L.values())
        APL = {str(L): int(round(-10.0*(sum_log10_P_L[L] - max_sum))) for L in cands}

        # Genotype enumeration
        pairs = all_genotype_pairs(cands)

        # Prior (Admixture); if extreme is specified, override
        prior = mix_prior_for_locus(freq_tab, lid, cands, popmix)
        if lid in extreme_map:
            prior = apply_extreme(prior, extreme_map[lid])

        like = {}
        post = {}
        prior_g = {}
        for (L1,L2) in pairs:
            # Likelihood
            s = 0.0
            for d in per_read_log10PL:
                pr = 0.5*(10**d[L1]) + 0.5*(10**d[L2])
                s += math.log10(max(pr, 1e-300))
            like[(L1,L2)] = s
            # Prior (HW)
            if L1 == L2:
                pg = prior[L1]*prior[L1]
            else:
                pg = 2.0*prior[L1]*prior[L2]
            prior_g[(L1,L2)] = pg
            post[(L1,L2)] = s + math.log10(max(pg, 1e-300))

        # Pick best/second best (like and post ranked independently)
        def top2(d):
            arr = sorted(d.items(), key=lambda kv: kv[1], reverse=True)
            best = arr[0]
            second = arr[1] if len(arr)>1 else (arr[0][0], float("-inf"))
            return best, second

        (best_like, sec_like) = top2(like)
        (best_post, sec_post) = top2(post)
        (gtL1, gtL2) = best_post[0]
        post_best, post2 = best_post[1], sec_post[1]
        like_best, like2 = best_like[1], sec_like[1]

        # Calculate Δ, if second best is -inf, treat as Δ=+inf -> GQ=99
        delta = post_best - post2 if math.isfinite(post2) else float('inf')
        GQ = _cap_gq_from_delta(delta)

        # Calculate log_prior_best (HW prior based on best genotype)
        def _fmt(x):
            return x if math.isfinite(x) else "NA"

        row = {
            "locus_id": lid,
            "motif": motif,
            "n_reads": n_reads,
            "candidates": ",".join(str(x) for x in cands),
            "GT": f"{gtL1}/{gtL2}",
            "GQ": GQ,
            "like_best_log10": like_best,
            "like2_log10": _fmt(like2),
            "post_best_log10": post_best,
            "post2_log10": _fmt(post2),
            "log_prior_best": math.log10(max(prior_g[(gtL1,gtL2)], 1e-300)),
            "prior_src": (f"{freq_type or 'none'}|" + (("mix:"+args.popmix) if args.popmix else "ALL_or_uniform")),
            "prior_L_json": json.dumps({str(L): prior[L] for L in cands}),
            "APL_json": json.dumps(APL, ensure_ascii=False),
        }
        out_rows.append(row)

    # Write TSV
    cols = ["locus_id","motif","n_reads","candidates","GT","GQ",
            "like_best_log10","like2_log10","post_best_log10","post2_log10",
            "log_prior_best","prior_src","prior_L_json","APL_json"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

if __name__ == "__main__":
    main()
