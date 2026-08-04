from pathlib import Path


def test_default_demo_uses_sw():
    root = Path(__file__).resolve().parents[1]
    for script in (root / "scripts/run_demo.ps1", root / "scripts/run_demo.sh"):
        text = script.read_text(encoding="utf-8")
        assert "--likelihood-backend phmm" not in text
        assert "--likelihood-model" not in text
