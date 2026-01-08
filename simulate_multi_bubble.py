import random
import os

# ================= 配置区域 =================
OUTPUT_PREFIX = "sim_data"
READ_LENGTH = 150       
INSERT_SIZE = 350       
INSERT_STD = 30         
SEQ_ERROR_RATE = 0.001  

# 【核心修正】
# 将侧翼长度增加到 1000bp。
# 之前 150bp 太短，导致 InsertSize(350) 的双端测序刚好“跨过”中间的 STR，
# 形成了 R1 在左、R2 在右，中间 STR 没测到的尴尬局面。
def generate_random_dna(length):
    return "".join(random.choice("ACGT") for _ in range(length))

random.seed(2024) # 固定种子
FLANK_L = generate_random_dna(1000)
FLANK_R = generate_random_dna(1000)

MOTIF = "CAG"
REF_REPEATS = 18 
CHROM_NAME = "chrSim"
BUILD_NAME = "SimRef_v1"

# 混合模型
SAMPLE_MIXTURE = [
    {
        "name": "Normal_Cells",
        "fraction": 0.60,       
        "alleles": [15, 18]  
    },
    {
        "name": "Tumor_Expanded",
        "fraction": 0.40,       
        "alleles": [15, 40]  
    }
]

# ================= 核心函数 =================
def generate_haplotype_sequence(n_repeats):
    return FLANK_L + (MOTIF * n_repeats) + FLANK_R

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def simulate_quality_scores(length):
    return "I" * length

def introduce_errors(seq):
    bases = ['A', 'C', 'G', 'T']
    seq_list = list(seq)
    for i in range(len(seq_list)):
        if random.random() < SEQ_ERROR_RATE:
            original = seq_list[i]
            choices = [b for b in bases if b != original]
            seq_list[i] = random.choice(choices)
    return "".join(seq_list)

def generate_reads(hap_seq, n_reads):
    r1_list, r2_list = [], []
    hap_len = len(hap_seq)
    
    # 为了保证覆盖度，我们多生成一点，因为现在基因组变大了
    # 但为了不改变用户参数，我们在循环内尝试随机剪切
    
    for _ in range(n_reads):
        # 随机插入片段长度
        current_insert = int(random.gauss(INSERT_SIZE, INSERT_STD))
        if current_insert > hap_len: current_insert = hap_len
        if current_insert < READ_LENGTH: current_insert = READ_LENGTH
            
        # 随机起始位置
        max_start = hap_len - current_insert
        if max_start < 0: max_start = 0
        start_pos = random.randint(0, max_start)
        end_pos = start_pos + current_insert
        fragment = hap_seq[start_pos:end_pos]
        
        # R1 (5'端)
        r1 = fragment[:READ_LENGTH]
        if len(r1) < READ_LENGTH: r1 = r1.ljust(READ_LENGTH, 'N')
            
        # R2 (3'端，反向互补)
        r2_raw = fragment[-READ_LENGTH:]
        if len(r2_raw) < READ_LENGTH: r2_raw = r2_raw.rjust(READ_LENGTH, 'N')
        r2 = reverse_complement(r2_raw)
        
        r1_list.append(introduce_errors(r1))
        r2_list.append(introduce_errors(r2))
    return r1_list, r2_list

def main():
    # 1. 生成参考基因组
    print(f"[*] 生成模拟参考序列 (Flank=1000bp)...")
    ref_seq = generate_haplotype_sequence(REF_REPEATS)
    ref_file = "sim_ref.fasta"
    with open(ref_file, "w") as f:
        f.write(f">{CHROM_NAME}\n{ref_seq}\n")
    print(f"   -> {ref_file}")

    # 2. 生成 STR Catalog
    start_pos = len(FLANK_L) + 1 
    end_pos = start_pos + (REF_REPEATS * len(MOTIF)) - 1
    
    catalog_file = "sim_catalog.tsv"
    with open(catalog_file, "w") as f:
        header = "chrom\tstart\tend\tmotif\tref_len\tlocus_id\tbuild"
        f.write(header + "\n")
        line = f"{CHROM_NAME}\t{start_pos}\t{end_pos}\t{MOTIF}\t{REF_REPEATS}\tHTT_Sim\t{BUILD_NAME}"
        f.write(line + "\n")
    print(f"   -> {catalog_file}")

    # 3. 生成 Reads
    # 增加 Reads 数量以适应更长的基因组 (5000 -> 10000)
    total_reads = 10000 
    print(f"[*] 生成混合样本 Reads (Total={total_reads})...")
    
    all_r1, all_r2 = [], []
    
    for comp in SAMPLE_MIXTURE:
        n = int(total_reads * comp["fraction"])
        print(f"   -> {comp['name']}: {n} reads, Alleles={comp['alleles']}")
        for repeats in comp["alleles"]:
            seq = generate_haplotype_sequence(repeats)
            r1, r2 = generate_reads(seq, n // 2)
            all_r1.extend(r1)
            all_r2.extend(r2)
            
    # Shuffle
    combined = list(zip(all_r1, all_r2))
    random.shuffle(combined)
    all_r1[:], all_r2[:] = zip(*combined)
    
    r1_name = f"{OUTPUT_PREFIX}_R1.fastq"
    r2_name = f"{OUTPUT_PREFIX}_R2.fastq"
    
    with open(r1_name, 'w') as f1, open(r2_name, 'w') as f2:
        for i, (s1, s2) in enumerate(zip(all_r1, all_r2)):
            h = f"@SIM_{i}"
            q = simulate_quality_scores(READ_LENGTH)
            f1.write(f"{h}/1\n{s1}\n+\n{q}\n")
            f2.write(f"{h}/2\n{s2}\n+\n{q}\n")
            
    print("[*] 模拟完成。")

if __name__ == "__main__":
    main()