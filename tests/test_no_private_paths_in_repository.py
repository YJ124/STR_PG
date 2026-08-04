from pathlib import Path


def test_no_private_paths_in_repository():
    root = Path(__file__).resolve().parents[1]
    forbidden = (
        "D:" + "\\STR-PG",
        "C:" + "\\Users\\Administrator",
        "sim-" + "data",
    )
    extensions = {".py", ".md", ".toml", ".yml", ".yaml", ".cff", ".txt", ".ps1", ".sh"}
    violations = []
    for path in root.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in extensions:
            continue
        relative = path.relative_to(root)
        # The formal local experiment workspace records external inputs,
        # interpreter paths, logs, and resumable commands. It is intentionally
        # excluded from the public source-release payload.
        if relative.parts[:1] == ("xiaorong",):
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        if any(value.lower() in text.lower() for value in forbidden):
            violations.append(str(relative))
    assert not violations
