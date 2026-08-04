$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $PSScriptRoot
Set-Location $ProjectRoot
New-Item -ItemType Directory -Path "demo_output" -Force | Out-Null

python -m strpg.cli build --ref examples/reference.fa --str-catalog examples/str_catalog.tsv --vcf examples/variants.vcf --out demo_output/graph.gfa --registry demo_output/registry.jsonl
python -m strpg.cli index --registry demo_output/registry.jsonl --out demo_output/index.sqlite
python -m strpg.cli map --index demo_output/index.sqlite --reads1 examples/demo_R1.fastq --reads2 examples/demo_R2.fastq --out demo_output/demo.gaf --tau0 0.05
cmd /d /c "python -m strpg.cli genotype --registry demo_output/registry.jsonl --gaf demo_output/demo.gaf --fq1 examples/demo_R1.fastq --fq2 examples/demo_R2.fastq --sample-id demo --read-cache-dir demo_output/sw_read_cache --out demo_output/genotypes.sw.tsv 2> demo_output/genotype.sw.stats.json"
if ($LASTEXITCODE -ne 0) { throw "SW genotype demo failed with exit code $LASTEXITCODE" }
