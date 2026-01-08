STR_PG: Pan-Genome Graph STR Genotyper
STR_PG is a comprehensive toolkit for constructing pan-genome graphs, mapping sequencing reads, and genotyping Short Tandem Repeats (STRs) using graph-based alignment and probabilistic modeling. It is designed to handle complex variation representation and provide accurate genotyping using population priors.

Features
Graph Construction: Build GFA-format graphs from reference genomes, VCF variants, and STR catalogs.

Efficient Indexing: Creates minimizer/syncmer seed indexes using DBM sharding for low-memory usage.

Fast Mapping: Multiprocessing-optimized read mapping to the graph (GAF output).

Probabilistic Genotyping: Genotype STRs using a Smith-Waterman alignment approach with population frequency priors.

Frequency Learning: Incrementally update population frequency databases based on genotyping results.

Requirements
System
OS: Linux / macOS

Python: 3.9+

Python Dependencies
Install the required packages using pip:

Bash

pip install pysam tqdm
pysam: Required for VCF parsing in the build step.

tqdm: Required for progress bars during mapping.

Standard libraries used: argparse, json, gzip, dbm, multiprocessing, pickle, struct.

Installation
Clone this repository:

Bash

git clone https://github.com/yourusername/STR_PG.git
cd STR_PG
Make the scripts executable (optional):

Bash

chmod +x *.py
Usage Pipeline
1. Build the Pan-Genome Graph (pgg_build.py)
Constructs a .gfa graph file containing the reference sequence, variants (SNP/Indel), and STR bubbles.

Input Requirements:

Reference Genome (FASTA)

STR Catalog (TSV)

Variants VCF (Optional)

STR Catalog Format (TSV): Columns: build, chrom, start, end, motif, candidate_Ls, pointer_id

Command:

Bash

python pgg_build.py \
    --ref reference.fa \
    --vcf variants.vcf.gz \
    --str-catalog str_catalog.tsv \
    --build GRCh38 \
    --flank 100 \
    --out graph.gfa
2. Index the Graph (pgg_index.py)
Creates a sharded index for fast read mapping. This generates path coordinates and seed databases.

Command:

Bash

python pgg_index.py \
    --graph graph.gfa \
    --out index_directory \
    --method syncmer \
    --k 15 --s 5 --t 2 \
    --shards 64
3. Map Reads to Graph (pgg_map_optimized.py)
Aligns FASTQ reads to the graph index. This step supports multi-processing for speed.

Command:

Bash

python pgg_map_optimized.py \
    --index index_directory \
    --reads1 sample_R1.fastq.gz \
    --reads2 sample_R2.fastq.gz \
    --out aligned_reads.gaf \
    --threads 8 \
    --batch_size 5000
Note: Use --reads for single-end sequencing.

4. Genotype STRs (pgg_genotype.py)
Performs the final genotyping using the Graph, Alignments (GAF), and raw reads.

Command:

Bash

python pgg_genotype.py \
    --gfa graph.gfa \
    --gaf aligned_reads.gaf \
    --fq1 sample_R1.fastq.gz \
    --fq2 sample_R2.fastq.gz \
    --freq freq_database.jsonl \
    --out genotypes.tsv \
    --region chr19:40000000-50000000 \
    --popmix "EUR=0.6,AFR=0.4"
--freq: Population frequency file (JSONL format).

--region: (Optional) Limit genotyping to a specific genomic region for speed.

--popmix: (Optional) Specify admixture proportions for priors.

5. (Optional) Update Frequency Database (pgg_update_freq.py)
Update your population frequency database using the genotyping results from a sample (uses homozygous loci only).

Command:

Bash

python pgg_update_freq.py \
    --geno genotypes.tsv \
    --freq-in current_freq.jsonl \
    --freq-out updated_freq.jsonl \
    --pop EAS \
    --min-GQ 20
6. (Optional) Validation/Testing (pgg_genotype_str_fix2.py)
A simplified genotyper useful for debugging specific loci or verifying priors without the full pipeline overhead.

Command:

Bash

python pgg_genotype_str_fix2.py \
    --gfa graph.gfa \
    --fq1 sample_R1.fastq.gz \
    --fq2 sample_R2.fastq.gz \
    --freq freq_database.jsonl \
    --out validation_output.tsv
File Formats
STR Catalog (TSV) Example
代码段

build	chrom	start	end	motif	candidate_Ls	pointer_id
GRCh38	chr1	1000	1004	AT	10,11,12,13	STR_1
Frequency Database (JSONL) Example
JSON

{"pointer_id": "STR_1", "freq": {"EAS": {"10": 0.2, "11": 0.8}, "EUR": {"10": 0.5, "11": 0.5}}}
License
MIT License (or specify your license here)
