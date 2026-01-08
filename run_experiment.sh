#!/bin/bash
# 遇到错误立即停止
set -e
TOOL_DIR="."
OUT_DIR="experiment_output"
REF="sim_ref.fasta"
CATALOG="sim_catalog.tsv"
R1="sim_data_R1.fastq"
R2="sim_data_R2.fastq"

echo "######################################################"
echo "   STR-PG 最终修正版测试：1000bp 侧翼 + 随机覆盖"
echo "######################################################"

if [ -d "$OUT_DIR" ]; then rm -rf "$OUT_DIR"; fi
mkdir -p $OUT_DIR

# Step 0
echo ""
echo "[Step 0] 生成模拟数据..."
python simulate_multi_bubble.py

# Step 1
echo ""
echo "[Step 1] Build..."
# 增加 flank 到 200，充分利用长侧翼
python $TOOL_DIR/pgg_build.py --ref $REF --str-catalog $CATALOG --build "SimRef_v1" --flank 200 --out $OUT_DIR/graph.gfa

# Step 2
echo ""
echo "[Step 2] Index..."
# Anchor 100 现在是安全的，因为我们的侧翼有 1000bp
python $TOOL_DIR/pgg_index.py --graph $OUT_DIR/graph.gfa --out $OUT_DIR/index --method syncmer --k 15 --w 10 --shards 16

# Step 3
echo ""
echo "[Step 3] Map..."
python $TOOL_DIR/pgg_map_optimized.py --index $OUT_DIR/index --reads1 $R1 --reads2 $R2 --out $OUT_DIR/mapped.gaf --method syncmer --threads 4 --max_occ 500 --slack 20

# Step 4
echo ""
echo "[Step 4] Genotype..."
# 增加 target_good_reads，因为现在reads总数多了
python $TOOL_DIR/pgg_genotype_optimized_bubble.py --gfa $OUT_DIR/graph.gfa --gaf $OUT_DIR/mapped.gaf --fq1 $R1 --fq2 $R2 --out $OUT_DIR/genotypes.tsv --target_good_reads 500 --max_reads_per_locus 5000

# Step 5
echo ""
echo "[Step 5] Verify..."
python verify_experiment.py

echo "######################################################"
echo "   实验完成！请检查:"
echo "   1. genotypes.tsv (n_reads 应该 > 0)"
echo "   2. output_verification.png"
echo "######################################################"