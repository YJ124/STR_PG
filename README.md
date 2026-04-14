# STR-PG: Decoupling Allelic Complexity on Pangenome Graph for Scalable STR Genotyping

**STR-PG** is an STR-specialized pangenome framework for building pointer-node-based graphs, mapping short reads with syncmer seeds, and genotyping short tandem repeats (STRs) using a **motif-aware profile hidden Markov model (pHMM)** and **population-aware Bayesian inference**.

Unlike conventional graph methods that explicitly encode every STR allele as a separate path, STR-PG represents each target STR locus as a **fixed-topology pointer node** in the graph and stores allele-resolved sequence content in an **external allele registry**. This design keeps local graph complexity stable, enables incremental allele registration without whole-graph reconstruction, and supports scalable short-read STR genotyping.

---

## Overview

STR-PG consists of four tightly connected components:

1. **Build**  
   Construct a GFA pangenome graph in which hypervariable STR loci are collapsed into **pointer nodes** rather than expanded into dense local bubbles.

2. **Index / Map**  
   Build a **syncmer-based graph index** and map reads to graph loci using seeding, chaining, and path-interval localization.

3. **Genotype**  
   Retrieve the candidate allele set from the external registry and perform probabilistic genotyping using:
   - population-aware priors
   - motif-aware pHMM likelihoods
   - posterior-based genotype selection

4. **Registry / Update**  
   Maintain a dynamic allele registry containing repeat counts, sequences, motif annotations, and population frequency priors, and update it incrementally without rebuilding the graph.

---

## Key Features

- **Pointer-node representation for STR loci**  
  Avoids topological “hairball” expansion caused by explicitly embedding many STR alleles in graph structure.

- **External allele registry**  
  Stores allele sequences, repeat counts, motif information, and population frequencies outside the graph.

- **Syncmer-based read localization**  
  Uses flank-derived syncmer seeds for coarse mapping, chaining, and locus assignment.

- **Motif-aware pHMM genotyping**  
  Uses a profile hidden Markov model with repeat-aware indel modeling and forward likelihood computation for STR allele evaluation.

- **Population-aware Bayesian inference**  
  Combines allele-frequency priors with per-read likelihoods to infer maximum a posteriori diploid genotypes.

- **GT / GQ output**  
  Reports genotype calls and posterior-derived genotype quality.

- **Incremental registry updates**  
  Supports metadata-level updates for newly supported alleles without whole-graph reconstruction.

- **Heterogeneity-aware design**  
  The framework is designed to preserve evidence for broadened or non-diploid local STR signals (“multi-bubbles”).

---

## Repository Contents

- `pgg_build.py`  
  Build STR-PG graph (`.gfa`) from reference sequence, STR catalog, and optional VCF.

- `pgg_index.py`  
  Build syncmer/minimizer-style graph index for fast mapping.

- `pgg_map_optimized.py`  
  Map single-end or paired-end reads to the graph and output GAF alignments.

- `pgg_genotype.py`  
  Main STR genotyper using motif-aware pHMM likelihoods and population priors.

- `pgg_update_freq.py`  
  Update population frequency / registry support information from high-confidence genotype calls.

- `pgg_genotype_hybrid_fast.py` *(optional, if included in your repo)*  
  A speed-oriented hybrid version that uses fast candidate screening before pHMM refinement.

- `pgg_genotype_str_fix2.py` *(optional / legacy / debug use)*  
  Lightweight validation script for targeted or debugging workflows.

---

## Requirements

### System
- Linux / macOS
- Python 3.9+

### Python dependencies

Install required packages:

```bash
pip install pysam tqdm

---

## Installation

Clone this repository:

```bash
git clone https://github.com/yourusername/STR_PG.git
cd STR_PG

```

Make the scripts executable (optional):

```bash
chmod +x *.py

```

---

## Usage Pipeline

### 1. Build the Pan-Genome Graph (`pgg_build.py`)

Constructs a `.gfa` graph file containing the reference sequence, variants (SNP/Indel), and STR bubbles.

**Input Requirements:**

* Reference Genome (FASTA)
* STR Catalog (TSV)
* Variants VCF (Optional)

**STR Catalog Format (TSV):**
Columns: `build`, `chrom`, `start`, `end`, `motif`, `candidate_Ls`, `pointer_id`

**Command:**

```bash
python pgg_build.py \
    --ref reference.fa \
    --vcf variants.vcf.gz \
    --str-catalog str_catalog.tsv \
    --build GRCh38 \
    --flank 100 \
    --out graph.gfa

```

The command to execute in the test_data folder is:
```bash
python pgg_build.py \
    --ref test_data/chr19.fa \
    --vcf test_data/chr19_STR.vcf.gz \
    --str-catalog test_data/str_catalog.tsv \
    --build GRCh38 \
    --flank 100 \
    --out test_data/graph.gfa

```
### 2. Index the Graph (`pgg_index.py`)

Creates a sharded index for fast read mapping. This generates path coordinates and seed databases.

**Command:**

```bash
python pgg_index.py \
    --graph test_data/graph.gfa \
    --out index_directory \
    --method syncmer \
    --k 15 --s 5 --t 2 \
    --shards 64

```

### 3. Map Reads to Graph (`pgg_map_optimized.py`)

Aligns FASTQ reads to the graph index. This step supports multi-processing for speed.

**Command:**

```bash
python pgg_map_optimized.py \
    --index index_directory \
    --reads1 sample_1_R1.fastq.gz \
    --reads2 sample_1_R2.fastq.gz \
    --out aligned_reads.gaf \
    --threads 8 \
    --batch_size 5000

```

The command to execute in the test_data folder is:
```bash
python pgg_map_optimized.py \
    --index index_dir \
    --reads1 test_data/sample_1_R1.fastq \
    --reads2 test_data/sample_1_R2.fastq \
    --out test_data/sample1.gaf \
    --threads 8 \
    --batch_size 5000

```
* *Note: Use `--reads` for single-end sequencing.*

### 4. Genotype STRs with motif-aware pHMM

localizes reads to pointer-node loci using graph alignments
retrieves candidate alleles from the registry
computes per-read likelihoods with a motif-aware pHMM
combines likelihoods with population priors
outputs MAP genotype (GT) and genotype quality (GQ)

**Command:**

```bash
python pgg_genotype.py \
    --gfa graph.gfa \
    --gaf aligned_reads.gaf \
    --fq1 sample_R1.fastq.gz \
    --fq2 sample_R2.fastq.gz \
    --freq freq_database.jsonl \
    --out genotypes.tsv \
    --region chr19:40000000-50000000 \
    --popmix "EUR=0.6,AFR=0.4"
```
Example in test_data
```bash
python pgg_genotype.py \
    --gfa test_data/graph.gfa \
    --gaf test_data/sample1.gaf \
    --fq1 test_data/sample_1_R1.fastq.gz \
    --fq2 test_data/sample_1_R2.fastq.gz \
    --freq freq19.jsonl \
    --out test_data/genotypes.tsv \
    --region chr19:40000000-50000000 \
    --popmix "EUR=0.6,AFR=0.4"

```
### 5. (Optional) Update Frequency Database (`pgg_update_freq.py`)

Update your population frequency database using the genotyping results from a sample (uses homozygous loci only).

**Command:**

```bash
python pgg_update_freq.py \
    --geno genotypes.tsv \
    --freq-in current_freq.jsonl \
    --freq-out updated_freq.jsonl \
    --pop EAS \
    --min-GQ 20

```
The command to execute in the test_data folder is:
```bash
python pgg_update_freq.py \
    --geno genotypes.tsv \
    --freq-in freq19.jsonl \
    --freq-out freq19_updated.jsonl \
    --pop EAS \
    --min-GQ 20

```
### 6. (Optional) Validation/Testing (`pgg_genotype_str_fix2.py`)

A simplified genotyper useful for debugging specific loci or verifying priors without the full pipeline overhead.

**Command:**

```bash
python pgg_genotype_str_fix2.py \
    --gfa graph.gfa \
    --fq1 sample_R1.fastq.gz \
    --fq2 sample_R2.fastq.gz \
    --freq freq_database.jsonl \
    --out validation_output.tsv

```

The command to execute in the test_data folder is:
```bash
python pgg_genotype_str_fix2.py \
    --gfa test_data/graph.gfa \
    --fq1 test_data/sample_1_R1.fastq.gz \
    --fq2 test_data/sample_1_R2.fastq.gz \
    --freq freq19.jsonl \
    --out validation_output.tsv

```
---

## File Formats

### STR Catalog (TSV) Example

```tsv
build	chrom	start	end	motif	candidate_Ls	pointer_id
GRCh38	chr1	1000	1004	AT	10,11,12,13	STR_1

```

### Frequency Database (JSONL) Example

```json
{"pointer_id": "STR_1", "freq": {"EAS": {"10": 0.2, "11": 0.8}, "EUR": {"10": 0.5, "11": 0.5}}}

```

