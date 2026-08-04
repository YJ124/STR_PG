from pathlib import Path
from strpg.graph import build_graph_and_registry
from strpg.index import build_seed_index
from strpg.mapper import map_fastq
from strpg.simulate import simulate_reads
from strpg.genotyper import genotype_dataset


def test_small_pipeline(tmp_path: Path):
    ref = tmp_path/"ref.fa"
    flank1=("ACGTCGATGCTAGCTACGAT"*15)[:250]
    flank2=("TGCACTGATCGTACGCTAGC"*15)[:250]
    seq=flank1+"CAG"*10+flank2
    ref.write_text(">chr1\n"+seq+"\n")
    cat=tmp_path/"catalog.tsv"
    cat.write_text(f"build\tchrom\tstart\tend\tmotif\tcandidate_Ls\tpointer_id\nGRCh38\tchr1\t{len(flank1)+1}\t{len(flank1)+30}\tCAG\t8,10,12\tP1\n")
    gfa,reg,idx=tmp_path/"g.gfa",tmp_path/"r.jsonl",tmp_path/"i.sqlite"
    build_graph_and_registry(ref,cat,gfa,reg,build="GRCh38",flank=100)
    build_seed_index(reg,idx)
    r1,r2=tmp_path/"R1.fq",tmp_path/"R2.fq"
    simulate_reads(reg,"P1","10=1",r1,r2,coverage=8,read_length=100,seed=5)
    gaf=tmp_path/"a.gaf"
    stats=map_fastq(idx,r1,gaf,r2,tau0=0.05)
    assert stats["mapped"]>0
    out=tmp_path/"gt.tsv"
    gstats=genotype_dataset(reg,gaf,r1,out,r2,"GLOBAL=1")
    assert gstats["calls"]==1
    assert "10/10" in out.read_text()
