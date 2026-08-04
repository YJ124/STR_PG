from pathlib import Path
from strpg.graph import build_graph_and_registry
from strpg.registry import AlleleRegistry


def test_pointer_graph_has_no_str_allele_paths(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\n" + "ACGT"*100 + "CAG"*5 + "TGCA"*100 + "\n")
    start = 401
    end = 415
    cat = tmp_path / "catalog.tsv"
    cat.write_text("build\tchrom\tstart\tend\tmotif\tcandidate_Ls\tpointer_id\nGRCh38\tchr1\t401\t415\tCAG\t4,5,6\tP1\n")
    gfa, reg = tmp_path / "graph.gfa", tmp_path / "registry.jsonl"
    build_graph_and_registry(ref, cat, gfa, reg, build="GRCh38", flank=40)
    text = gfa.read_text()
    assert "TP:Z:STR_POINTER" in text
    assert "allele_" not in text
    assert "L4" not in text and "L5" not in text and "L6" not in text
    registry = AlleleRegistry.load(reg)
    assert [a.repeat_count for a in registry.loci["P1"].alleles] == [4,5,6]
