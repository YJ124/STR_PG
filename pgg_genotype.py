#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pgg_genotype.py
---------------------------------------------------------
Author:YJ

Includes:
1. Paired-end Joint Calling
2. Region Filter - Core speed-up feature
3. Performance parameter optimization (Band=50, MaxReads=100)
4. Complete frequency loading and prior calculation logic

Usage:
python pgg_genotype_v13.py --gfa ... --gaf ... --fq1 ... --fq2 ... --out ... \
    --region chr19:40000000-50000000
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

# -------------------------------------------------------
# Global Settings and Tools
# -------------------------------------------------------
random.seed(1)
TRANS_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def revcomp(seq):
    return seq.translate(TRANS_TABLE)[::-1]

# -------------------------------------------------------
# ID Cleaning and Coordinate Mapping
# -------------------------------------------------------
def normalize_read_id(raw_id):
    """
    Clean Read ID, remove descriptions after /1, /2 or spaces, keep core ID for pairing.
    """
    if raw_id.startswith("@"): 
        raw_id = raw_id[1:]
    core = raw_id.split()[0] # Take the part before space
    # Remove trailing pairing indicators
    if core.endswith("/1") or core.endswith("/2"):
        core = core[:-2]
    elif core.endswith(".1") or core.endswith(".2"):
        core = core[:-2]
    return core

def parse_region(region_str):
    """
    Parse region argument, format chr:start-end
    """
    if not region_str: return None
    try:
        chrom, ran = region_str.split(":")
        start, end = map(int, ran.split("-"))
        return (chrom, start, end)
    except:
        return None

def is_in_region(locus_id, region_tuple):
    """
    Check if locus_id is within the specified region
    locus_id format example: GRCh38:chr19:3207858-3207897:CTG
    """
    if not region_tuple: return True
    try:
        parts = locus_id.split(":")
        # Find parts containing coordinates
        l_chrom, l_start, l_end = None, 0, 0
        
        # Robust parsing: Look for fields containing 'chr' and numeric fields containing '-'
        for p in parts:
            if "chr" in p: 
                l_chrom = p
            elif "-" in p and p.replace("-", "").isdigit():
                s, e = map(int, p.split("-"))
                l_start, l_end = s, e
        
        # If parsing fails, return False conservatively to avoid deviations
        if l_chrom is None or l_start == 0:
            return False

        r_chrom, r_start, r_end = region_tuple
        
        # Chromosomes must match
        if l_chrom != r_chrom: return False
        
        # Check for overlap
        # locus: [l_start, l_end], region: [r_start, r_end]
        if l_end < r_start or l_start > r_end:
            return False
            
        return True
    except:
        return False

def parse_coords_from_gfa_id(locus_id):
    try:
        parts = locus_id.split(':')
        for i, p in enumerate(parts):
            if '-' in p and p.replace('-', '').isdigit():
                return f"{parts[i-1]}:{p}"
    except: pass
    return None

def parse_coords_from_gaf_id(path_name):
    try:
        if path_name.startswith("allele_"):
            parts = path_name.split('_')
            # Assume path_name format allele_chr19_123_456_...
            # Needs adjustment based on actual naming rules during GFA build
            # Attempting a generalization here
            if len(parts) >= 4:
                return f"{parts[1]}:{parts[2]}-{parts[3]}"
    except: pass
    return None

# -------------------------------------------------------
# Helper Functions: Probability and GQ
# -------------------------------------------------------
def _cap_gq_from_delta(delta):
    if not math.isfinite(delta): return 99
    if delta < 0: delta = 0.0
    delta = min(delta, 9.9)
    return int(round(10.0 * delta))

def parse_popmix(s):
    if not s: return {}
    out = {}
    for kv in s.split(","):
        if not kv.strip(): continue
        try:
            k, v = kv.split("=")
            out[k.strip().upper()] = float(v)
        except: pass
    tot = sum(out.values())
    if tot > 0: 
        for k in out: 
            out[k] /= tot
    return out

def parse_extreme_list(lst):
    out = {}
    for s in lst or []:
        if "=" in s:
            lid, Ls = s.split("=", 1)
            try: out[lid.strip()] = int(Ls)
            except: pass
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
            for p, w in popmix.items(): s += w * pops.get(p, 0.0)
            raw[L] = s
        else:
            raw[L] = pops.get("ALL", 0.0)
        if raw[L] > 0: has_any = True
    
    if not has_any:
        for L in cand_L: 
            raw[L] = 1.0
            
    tot = sum(raw.values()) + len(cand_L) * eps
    return {L: (v + eps) / tot for L, v in raw.items()}

def apply_extreme(prior, limit_val):
    # Simple implementation: Set the probability of extreme L to be very high
    new_prior = {}
    eps = 1e-6
    n = len(prior)
    for L in prior:
        if L == limit_val:
            new_prior[L] = 1.0 - (n-1)*eps
        else:
            new_prior[L] = eps
    return new_prior

# -------------------------------------------------------
# I/O: GAF Parsing and FASTQ Extraction
# -------------------------------------------------------
def parse_gaf_and_build_index(gaf_path, loci_map):
    print(f"[INFO] Loading GAF index from: {gaf_path}", file=sys.stderr)
    coord_to_gfa_id = {}
    for lid in loci_map:
        coord_key = parse_coords_from_gfa_id(lid)
        if coord_key: coord_to_gfa_id[coord_key] = lid
            
    mapping = defaultdict(set)
    op = gzip.open if str(gaf_path).endswith(".gz") else open
    try:
        with op(gaf_path, "rt", encoding="utf-8", errors="ignore") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 6: continue
                raw_name = parts[0]
                path_name = parts[5]
                if path_name == "*" or not path_name: continue
                
                target_lid = None
                if path_name in loci_map:
                    target_lid = path_name
                else:
                    # Try matching via coordinates
                    c_key = parse_coords_from_gaf_id(path_name)
                    if c_key and c_key in coord_to_gfa_id:
                        target_lid = coord_to_gfa_id[c_key]
                
                if target_lid:
                    clean_id = normalize_read_id(raw_name)
                    mapping[clean_id].add(target_lid)
    except Exception as e:
        sys.exit(f"[ERROR] GAF Parsing Failed: {e}")
    print(f"[INFO] GAF Index ready. Found reads for {len(mapping)} unique IDs.", file=sys.stderr)
    return mapping

def stream_fastq_and_group(fq1, fq2, read_to_loci):
    """
    Read paired-end FASTQ, store R1 and R2 bound together.
    Return structure: Dict[locus_id] -> List of ( (n1, s1, q1), (n2, s2, q2) )
    """
    loci_data = defaultdict(list)
    op1 = gzip.open if str(fq1).endswith(".gz") else open
    op2 = gzip.open if str(fq2).endswith(".gz") else open
    print(f"[INFO] Streaming FASTQ pairs...", file=sys.stderr)
    
    total_pairs = 0
    kept_pairs = 0
    
    try:
        with op1(fq1, "rt") as f1, op2(fq2, "rt") as f2:
            while True:
                h1 = f1.readline().strip()
                if not h1: break # EOF
                s1 = f1.readline().strip()
                _  = f1.readline()
                q1 = f1.readline().strip()
                
                h2 = f2.readline().strip()
                s2 = f2.readline().strip()
                _  = f2.readline()
                q2 = f2.readline().strip()
                
                if not h2: break 

                total_pairs += 1
                if total_pairs % 5000000 == 0:
                    print(f"  Processed {total_pairs} pairs... (Kept {kept_pairs})", file=sys.stderr)
                
                core_id = normalize_read_id(h1)
                
                if core_id in read_to_loci:
                    kept_pairs += 1
                    target_loci = read_to_loci[core_id]
                    pair_data = ( (h1, s1, q1), (h2, s2, q2) )
                    for lid in target_loci:
                        loci_data[lid].append(pair_data)
                        
    except Exception as e: 
        sys.exit(f"[ERROR] FASTQ Reading Failed: {e}")
        
    print(f"[INFO] Extraction done. Kept {kept_pairs} relevant pairs.", file=sys.stderr)
    return loci_data

# -------------------------------------------------------
# Smith-Waterman & Likelihood Calculation
# -------------------------------------------------------
def sw_score(a, b, match=2, mismatch=-5, gap=-5, band=None): 
    n, m = len(a), len(b)
    prev = [0]*(m+1)
    best = 0
    
    for i in range(1, n+1):
        cur = [0]*(m+1)
        ai = a[i-1]
        
        j_lo, j_hi = 1, m
        if band:
            j_lo = max(1, i - band)
            j_hi = min(m, i + band)
            
        for j in range(j_lo, j_hi+1):
            sc = max(0, 
                     prev[j-1] + (match if ai == b[j-1] else mismatch),
                     prev[j] + gap, 
                     cur[j-1] + gap)
            cur[j] = sc
            if sc > best: best = sc
        prev = cur
    return best

def get_read_log10_likelihood(seq, H_dict, args):
    if not seq: return None
    
    qseq = seq
    qseq_rc = revcomp(seq)
    
    scores = {}
    bestS = -1.0
    
    # Optimization: Reduce Band
    for L, hseq in H_dict.items():
        s1 = sw_score(qseq, hseq, band=args.band)
        s2 = sw_score(qseq_rc, hseq, band=args.band)
        s = max(s1, s2)
        scores[L] = s
        if s > bestS: bestS = s
    
    # Filter low-quality alignments
    ideal_score = len(qseq) * 2
    if bestS < (ideal_score * args.min_score_ratio):
        return None
    
    log10P = {}
    denom = 0.0
    tau = max(args.tau, 1e-6)
    
    exp_vals = {}
    for L, s in scores.items():
        val = math.exp((s - bestS) / tau)
        exp_vals[L] = val
        denom += val
        
    for L in scores.keys():
        p = max(exp_vals[L] / denom, 1e-300)
        log10P[L] = math.log10(p)
        
    return log10P

# -------------------------------------------------------
# Configuration Parsing
# -------------------------------------------------------
def parse_gfa_strs(path):
    loci = []
    try:
        with open(path, "r") as f:
            for ln in f:
                if not ln or ln.startswith("#"): continue
                parts = ln.rstrip("\n").split("\t")
                if len(parts) >= 7 and parts[0] == "X" and parts[1] == "STR":
                    try: cand_L = [int(x) for x in parts[6].split(",") if x]
                    except: cand_L = []
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
    if not path: return tab, "none"
    if not os.path.exists(path):
        print(f"[WARN] Frequency file not found: {path}, using uniform prior.", file=sys.stderr)
        return tab, "missing"
    # Compatible with JSONL and JSON
    try:
        with open(path, 'r') as f:
            for line in f:
                if not line.strip(): continue
                try:
                    obj = json.loads(line)
                    # Support list of dicts (full JSON) or single dict (JSONL)
                    if isinstance(obj, list):
                        items = obj
                    else:
                        items = [obj]
                        
                    for item in items:
                        lid = item.get("locus_id") or item.get("id") or item.get("pointer_id")
                        if not lid: continue
                        data = item.get("L_freq") or item.get("freq") or item
                        for k, v in data.items():
                            if k.isdigit() and isinstance(v, dict):
                                tab[lid][int(k)] = v
                            elif k in ["AFR", "EUR", "EAS", "AMR", "SAS", "ALL"] and isinstance(v, dict):
                                # Reverse parsing is also possible
                                pass
                except: pass
        return tab, "json"
    except: return tab, "error"

# -------------------------------------------------------
# Main Logic
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
    ap.add_argument("--tau", type=float, default=12.0)
    
    # Optimization parameters
    ap.add_argument("--band", type=int, default=50) # Reduced to 50
    ap.add_argument("--target_good_reads", type=int, default=50) # Reduced to 50
    ap.add_argument("--max_reads_per_locus", type=int, default=100) # Reduced to 100
    ap.add_argument("--min_score_ratio", type=float, default=0.3)
    
    # [Core Parameter] Region Filter
    ap.add_argument("--region", required=False, help="Filter loci by region (e.g. chr19:40000000-50000000)")
    
    args = ap.parse_args()
    
    out_dir = os.path.dirname(args.out)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    popmix = parse_popmix(args.popmix)
    extreme_map = parse_extreme_list(args.extreme)

    # 1. Parse Region
    region_filter = parse_region(args.region)
    if region_filter:
        print(f"[INFO] Region Filter Enabled: {args.region}", file=sys.stderr)

    # 2. Parse GFA and apply filtering
    all_loci = parse_gfa_strs(args.gfa)
    loci_map = {}
    for l in all_loci:
        if is_in_region(l["locus_id"], region_filter):
            loci_map[l["locus_id"]] = l
    
    print(f"[INFO] Loci active after filtering: {len(loci_map)} (Original GFA: {len(all_loci)})", file=sys.stderr)
    if len(loci_map) == 0:
        sys.exit("[ERROR] No loci found in the specified region! Check your --region or GFA.")

    freq_tab, ftype = load_freq_auto(args.freq)
    
    # 3. Parse GAF to get whitelist (Only parse reads corresponding to whitelist Loci)
    read_whitelist = parse_gaf_and_build_index(args.gaf, loci_map)
    
    # 4. Extract and group FASTQ
    grouped_reads = stream_fastq_and_group(args.fq1, args.fq2, read_whitelist)

    out_rows = []
    
    total_active_loci = len(grouped_reads)
    print(f"[INFO] Genotyping {total_active_loci} active loci with reads...", file=sys.stderr)

    start_time = time.time()
    processed_count = 0

    # 5. Genotype locus by locus
    for lid in loci_map:
        if lid not in grouped_reads:
            # If empty results are desired, that's fine too; here we only output those with reads for brevity
            # out_rows.append({"locus_id": lid, "GT": "./.", "GQ": 0, "n_reads": 0})
            continue
        
        processed_count += 1
        if processed_count % 100 == 0:
            elapsed = time.time() - start_time
            print(f"  [Progress] {processed_count}/{total_active_loci} loci processed ({elapsed:.1f}s)...", file=sys.stderr)

        pairs_for_this_locus = grouped_reads[lid]
        loc = loci_map[lid]
        cands = sorted(set(loc["candidates"])) 
        if not cands: cands = list(range(3, 20))
        
        H = {L: (loc["left"] + loc["motif"]*L + loc["right"]) for L in cands}
        
        per_pair_joint_log10P = []
        valid_pairs_count = 0
        total_checked = 0
        
        for (r1_data, r2_data) in pairs_for_this_locus:
            if valid_pairs_count >= args.target_good_reads: break
            if total_checked >= args.max_reads_per_locus: break
            total_checked += 1

            (n1, s1, q1) = r1_data
            (n2, s2, q2) = r2_data
            
            logP_R1 = get_read_log10_likelihood(s1, H, args)
            logP_R2 = get_read_log10_likelihood(s2, H, args)
            
            if logP_R1 is None and logP_R2 is None: continue 
            
            valid_pairs_count += 1
            
            joint_logP = {}
            for L in cands:
                s = 0.0
                if logP_R1: s += logP_R1[L]
                if logP_R2: s += logP_R2[L]
                joint_logP[L] = s
            
            per_pair_joint_log10P.append(joint_logP)

        if valid_pairs_count == 0:
            out_rows.append({"locus_id": lid, "GT": "./.", "GQ": 0, "n_reads": 0})
            continue

        prior = mix_prior_for_locus(freq_tab, lid, cands, popmix)
        if lid in extreme_map: prior = apply_extreme(prior, extreme_map[lid])
        
        post = {}
        cand_pairs = []
        for i, L1 in enumerate(cands):
            for L2 in cands[i:]:
                cand_pairs.append((L1,L2))
                
        for (L1,L2) in cand_pairs:
            ll = 0.0
            for d in per_pair_joint_log10P:
                p_l1 = 10**d[L1]
                p_l2 = 10**d[L2]
                pr = 0.5 * p_l1 + 0.5 * p_l2
                ll += math.log10(max(pr, 1e-300))
            
            pg = (prior[L1]**2) if L1==L2 else (2*prior[L1]*prior[L2])
            post[(L1,L2)] = ll + math.log10(max(pg, 1e-300))
            
        sorted_post = sorted(post.items(), key=lambda x: x[1], reverse=True)
        (gtL1, gtL2), best_p = sorted_post[0]
        second_p = sorted_post[1][1] if len(sorted_post)>1 else -999.0
        
        GQ = _cap_gq_from_delta(best_p - second_p)
        
        out_rows.append({
            "locus_id": lid, 
            "GT": f"{gtL1}/{gtL2}", 
            "GQ": GQ, 
            "n_reads": valid_pairs_count
        })

    # Output results
    keys = ["locus_id","GT","GQ","n_reads"]
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
