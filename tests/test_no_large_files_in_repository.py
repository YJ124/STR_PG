from pathlib import Path


def test_no_large_files_in_repository():
    root = Path(__file__).resolve().parents[1]
    ignored = {".git", ".venv", "demo_output", "read_cache", "xiaorong"}
    large = [
        path.relative_to(root)
        for path in root.rglob("*")
        if path.is_file()
        and not ignored.intersection(path.relative_to(root).parts)
        and path.stat().st_size > 2 * 1024 * 1024
    ]
    assert not large
