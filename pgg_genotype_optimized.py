#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pgg_genotype_v12_final.py
---------------------------------------------------------
【最终完整版 V12】
结合了 V10 的正确双端逻辑 + V11 的性能安全锁。

功能特性：
1. 双端联合检测 (Joint Calling): 修复了准确率倒挂问题。
2. 性能保护 (Safety Brake): 防止 GAF 多重比对导致的计算爆炸。
3. 实时监控: 每处理 100 个位点打印一次进度。

使用方法：
python pgg_genotype_v12_final.py --gfa ... --gaf ... --fq1 ... --fq2 ... --out ... \
    --target_good_reads 200 --max_reads_per_locus 1000
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
# 全局设置与工具
# -------------------------------------------------------
random.seed(1)
TRANS_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def revcomp(seq):
    return seq.translate(TRANS_TABLE)[::-1]

# -------------------------------------------------------
# ID 清洗与坐标映射
# -------------------------------------------------------
def normalize_read_id(raw_id):
    """
    清洗 Read ID，去除 /1, /2 或空格后的描述，保留核心 ID 以便配对。
    """
    if raw_id.startswith("@"): 
        raw_id = raw_id[1:]
    core = raw_id.split()[0] # 取空格前部分
    # 去除末尾的配对标识
    if core.endswith("/1") or core.endswith("/2"):
        core = core[:-2]
    elif core.endswith(".1") or core.endswith(".2"):
        core = core[:-2]
    return core

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
            if len(parts) >= 4:
                return f"{parts[1]}:{parts[2]}-{parts[3]}"
    except: pass
    return None

# -------------------------------------------------------
# 辅助函数：概率与 GQ
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
    return prior

# -------------------------------------------------------
# I/O：GAF 解析与 FASTQ 提取 (核心修复)
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
                    c_key = parse_coords_from_gaf_id(path_name)
                    if c_key and c_key in coord_to_gfa_id:
                        target_lid = coord_to_gfa_id[c_key]
                
                if target_lid:
                    # GAF 中的 read name 可能带也可能不带 /1 /2，统一清洗
                    clean_id = normalize_read_id(raw_name)
                    mapping[clean_id].add(target_lid)
    except Exception as e:
        sys.exit(f"[ERROR] GAF Parsing Failed: {e}")
    print(f"[INFO] GAF Index ready. Found reads for {len(mapping)} unique IDs.", file=sys.stderr)
    return mapping

def stream_fastq_and_group(fq1, fq2, read_to_loci):
    """
    【核心修复】
    读取双端 FASTQ，将 R1 和 R2 绑定在一起存储，而不是拆散。
    返回结构: Dict[locus_id] -> List of ( (n1, s1, q1), (n2, s2, q2) )
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
                # 同时读取 R1 和 R2 的 4 行
                h1 = f1.readline().strip()
                if not h1: break # EOF
                s1 = f1.readline().strip()
                _  = f1.readline()
                q1 = f1.readline().strip()
                
                h2 = f2.readline().strip()
                s2 = f2.readline().strip()
                _  = f2.readline()
                q2 = f2.readline().strip()
                
                # 确保文件长度一致
                if not h2: break 

                total_pairs += 1
                if total_pairs % 2000000 == 0:
                    print(f"  Processed {total_pairs} pairs... (Kept {kept_pairs})", file=sys.stderr)
                
                # 使用清洗后的 Core ID 进行匹配
                core_id = normalize_read_id(h1)
                
                if core_id in read_to_loci:
                    kept_pairs += 1
                    target_loci = read_to_loci[core_id]
                    
                    # 打包成 Tuple，不拆散！
                    pair_data = ( (h1, s1, q1), (h2, s2, q2) )
                    
                    for lid in target_loci:
                        loci_data[lid].append(pair_data)
                        
    except Exception as e: 
        sys.exit(f"[ERROR] FASTQ Reading Failed: {e}")
        
    print(f"[INFO] Extraction done. Kept {kept_pairs} relevant pairs.", file=sys.stderr)
    return loci_data

# -------------------------------------------------------
# Smith-Waterman & 似然计算 (模块化)
# -------------------------------------------------------
def sw_score(a, b, match=2, mismatch=-5, gap=-5, band=None): 
    n, m = len(a), len(b)
    # 内存优化：只保留两行
    prev = [0]*(m+1)
    best = 0
    
    # 简单的 Banding 策略
    for i in range(1, n+1):
        cur = [0]*(m+1)
        ai = a[i-1]
        
        j_lo, j_hi = 1, m
        if band:
            j_lo = max(1, i - band)
            j_hi = min(m, i + band)
            
        for j in range(j_lo, j_hi+1):
            # 局部比对，得分 >= 0
            sc = max(0, 
                     prev[j-1] + (match if ai == b[j-1] else mismatch),
                     prev[j] + gap, 
                     cur[j-1] + gap)
            cur[j] = sc
            if sc > best: best = sc
        prev = cur
    return best

def get_read_log10_likelihood(seq, H_dict, args):
    """
    计算单条 Read 针对所有候选等位基因的 Log10 似然。
    如果 Read 质量太差（得分低），返回 None。
    """
    if not seq: return None
    
    qseq = seq
    qseq_rc = revcomp(seq)
    
    scores = {}
    bestS = -1.0
    
    # 1. SW 比对
    for L, hseq in H_dict.items():
        s1 = sw_score(qseq, hseq, band=args.band)
        s2 = sw_score(qseq_rc, hseq, band=args.band)
        s = max(s1, s2)
        scores[L] = s
        if s > bestS: bestS = s
    
    # 2. 垃圾过滤 (Dynamic Threshold)
    # 理想得分大约是 2 * read_len (全匹配)
    ideal_score = len(qseq) * 2
    if bestS < (ideal_score * args.min_score_ratio):
        return None # 标记为无效 Read
    
    # 3. Softmax 转概率
    log10P = {}
    denom = 0.0
    tau = max(args.tau, 1e-6)
    
    exp_vals = {}
    for L, s in scores.items():
        val = math.exp((s - bestS) / tau)
        exp_vals[L] = val
        denom += val
        
    # 4. 归一化并取对数
    for L in scores.keys():
        p = max(exp_vals[L] / denom, 1e-300)
        log10P[L] = math.log10(p)
        
    return log10P

# -------------------------------------------------------
# 主逻辑
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
    try:
        with open(path, 'r') as f:
            for line in f:
                if not line.strip(): continue
                try:
                    obj = json.loads(line)
                    lid = obj.get("locus_id") or obj.get("id")
                    if not lid: continue
                    data = obj.get("L_freq") or obj.get("freq") or obj
                    for k, v in data.items():
                        if k.isdigit() and isinstance(v, dict):
                            tab[lid][int(k)] = v
                except: pass
        return tab, "json"
    except: return tab, "error"

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
    ap.add_argument("--band", type=int, default=500) 
    
    # 性能控制参数
    ap.add_argument("--target_good_reads", type=int, default=200, help="Stop after finding N good pairs")
    ap.add_argument("--min_score_ratio", type=float, default=0.3, help="Filter reads with score < ratio * ideal") 
    
    # 【新增功能】：安全刹车参数
    ap.add_argument("--max_reads_per_locus", type=int, default=1000, 
                    help="Hard Limit: Stop checking a locus after processing this many pairs, even if target not met.")
    
    args = ap.parse_args()
    
    # 创建输出目录（如果不存在）
    out_dir = os.path.dirname(args.out)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    popmix = parse_popmix(args.popmix)
    extreme_map = parse_extreme_list(args.extreme)

    loci = parse_gfa_strs(args.gfa)
    if not loci: sys.exit("[ERROR] No STR loci found in GFA.")
    loci_map = {l["locus_id"]: l for l in loci}
    
    freq_tab, ftype = load_freq_auto(args.freq)
    
    # 1. 解析 GAF 获取白名单
    read_whitelist = parse_gaf_and_build_index(args.gaf, loci_map)
    
    # 2. 提取并分组 FASTQ (Key: locus_id, Value: List of Pairs)
    grouped_reads = stream_fastq_and_group(args.fq1, args.fq2, read_whitelist)

    out_rows = []
    
    total_active_loci = len(grouped_reads)
    print(f"[INFO] Genotyping {total_active_loci} active loci...", file=sys.stderr)

    # 计时开始
    start_time = time.time()
    processed_count = 0

    # 3. 逐个位点进行基因分型
    for lid in loci_map:
        # 如果这个位点没有任何 Reads，直接跳过
        if lid not in grouped_reads:
            out_rows.append({"locus_id": lid, "GT": "./.", "GQ": 0, "n_reads": 0})
            continue
        
        # 【进度条】：每处理 100 个位点打印一次，确保你知道程序还活着
        processed_count += 1
        if processed_count % 100 == 0:
            elapsed = time.time() - start_time
            print(f"  [Progress] {processed_count}/{total_active_loci} loci processed ({elapsed:.1f}s)...", file=sys.stderr)

        pairs_for_this_locus = grouped_reads[lid]
        loc = loci_map[lid]
        cands = sorted(set(loc["candidates"])) 
        if not cands: cands = list(range(3, 20))
        
        # 构建参考单倍型序列
        H = {L: (loc["left"] + loc["motif"]*L + loc["right"]) for L in cands}
        
        per_pair_joint_log10P = []
        valid_pairs_count = 0
        total_checked_pairs = 0 # 记录该位点总共检查了多少对
        
        # 遍历所有 Read Pairs
        for (r1_data, r2_data) in pairs_for_this_locus:
            # 停止条件1：找到足够的有效 Reads (质量达标)
            if valid_pairs_count >= args.target_good_reads:
                break
            
            # 【安全刹车】：停止条件2：硬性止损
            # 如果查了 max_reads_per_locus 条 reads 还没凑齐，强制停止，防止死循环计算
            if total_checked_pairs >= args.max_reads_per_locus:
                break
                
            total_checked_pairs += 1

            (n1, s1, q1) = r1_data
            (n2, s2, q2) = r2_data
            
            # 分别计算 R1 和 R2 的似然
            logP_R1 = get_read_log10_likelihood(s1, H, args)
            logP_R2 = get_read_log10_likelihood(s2, H, args)
            
            # 如果两个都是垃圾（None），直接跳过
            if logP_R1 is None and logP_R2 is None:
                continue 
            
            valid_pairs_count += 1
            
            joint_logP = {}
            for L in cands:
                # 联合似然计算：叠加 LogP
                score = 0.0
                if logP_R1 is not None and logP_R2 is not None:
                    score = logP_R1[L] + logP_R2[L]
                elif logP_R1 is not None:
                    score = logP_R1[L]
                elif logP_R2 is not None:
                    score = logP_R2[L]
                
                joint_logP[L] = score
            
            per_pair_joint_log10P.append(joint_logP)

        # 所有的 Pairs 处理完毕
        if valid_pairs_count == 0:
            out_rows.append({"locus_id": lid, "GT": "./.", "GQ": 0, "n_reads": 0})
            continue

        # 4. 贝叶斯基因分型 (Bayesian Genotyping)
        prior = mix_prior_for_locus(freq_tab, lid, cands, popmix)
        if lid in extreme_map: prior = apply_extreme(prior, extreme_map[lid])
        
        post = {}
        # 遍历所有可能的基因型组合 (L1, L2)
        cand_pairs = []
        for i, L1 in enumerate(cands):
            for L2 in cands[i:]:
                cand_pairs.append((L1,L2))
                
        for (L1,L2) in cand_pairs:
            # Log Likelihood = Sum Log(P(ReadPair_i | GT))
            ll = 0.0
            for d in per_pair_joint_log10P:
                log_p_l1 = d[L1]
                log_p_l2 = d[L2]
                
                # P(Data | L1,L2) = 0.5 * P(Data | L1) + 0.5 * P(Data | L2)
                p_l1 = 10**log_p_l1
                p_l2 = 10**log_p_l2
                pr = 0.5 * p_l1 + 0.5 * p_l2
                
                ll += math.log10(max(pr, 1e-300))
            
            # 乘以先验
            pg = (prior[L1]**2) if L1==L2 else (2*prior[L1]*prior[L2])
            post[(L1,L2)] = ll + math.log10(max(pg, 1e-300))
            
        # 5. 排序找最大后验
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

    # 6. 输出结果
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