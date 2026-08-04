from __future__ import annotations

import json
import sqlite3
from pathlib import Path

from .registry import AlleleRegistry
from .syncmer import iter_syncmers

SCHEMA = """
PRAGMA journal_mode=WAL;
PRAGMA synchronous=NORMAL;
CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS loci (
  pointer_id TEXT PRIMARY KEY,
  chrom TEXT NOT NULL,
  start INTEGER NOT NULL,
  end INTEGER NOT NULL,
  motif TEXT NOT NULL,
  anchor_length INTEGER NOT NULL,
  path_length INTEGER NOT NULL,
  flank_seed_count INTEGER NOT NULL,
  left_seed_count INTEGER NOT NULL,
  right_seed_count INTEGER NOT NULL
);
CREATE TABLE IF NOT EXISTS seeds (
  seed TEXT NOT NULL,
  pointer_id TEXT NOT NULL,
  side TEXT NOT NULL,
  path_pos INTEGER NOT NULL,
  orientation INTEGER NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_seed ON seeds(seed);
CREATE INDEX IF NOT EXISTS idx_seed_pointer ON seeds(seed, pointer_id);
"""


class SeedIndex:
    def __init__(self, path: str | Path, read_only: bool = True):
        uri = f"file:{Path(path).resolve()}?mode=ro" if read_only else str(path)
        self.conn = sqlite3.connect(uri, uri=read_only)
        self.conn.row_factory = sqlite3.Row
        self.metadata = {row["key"]: json.loads(row["value"]) for row in self.conn.execute("SELECT key,value FROM metadata")}
        # STR-PG indexes contain only locus-flank seeds.  Loading these rows once
        # avoids one SQLite statement for every read syncmer while preserving
        # the exact lookup and maximum-occurrence semantics.
        self._seed_cache: dict[str, list[sqlite3.Row]] = {}
        for row in self.conn.execute(
            "SELECT seed,pointer_id,side,path_pos,orientation FROM seeds"
        ):
            self._seed_cache.setdefault(row["seed"], []).append(row)
        self._locus_cache = {
            row["pointer_id"]: row for row in self.conn.execute("SELECT * FROM loci")
        }

    def close(self) -> None:
        self.conn.close()

    def lookup(self, seed: str, max_occurrences: int = 200) -> list[sqlite3.Row]:
        rows = self._seed_cache.get(seed, [])
        return [] if len(rows) > max_occurrences else rows

    def locus(self, pointer_id: str) -> sqlite3.Row:
        row = self._locus_cache.get(pointer_id)
        if row is None:
            raise KeyError(pointer_id)
        return row


def build_seed_index(
    registry_path: str | Path,
    out_path: str | Path,
    k: int = 15,
    s: int = 5,
    t: int = 2,
) -> None:
    out = Path(out_path)
    if out.exists():
        out.unlink()
    registry = AlleleRegistry.load(registry_path)
    conn = sqlite3.connect(out)
    conn.executescript(SCHEMA)
    metadata = {"version": 1, "method": "open_syncmer", "k": k, "s": s, "t": t, "coordinate_model": "left_flank + pointer(1bp) + right_flank"}
    conn.executemany("INSERT INTO metadata(key,value) VALUES(?,?)", [(key, json.dumps(value)) for key, value in metadata.items()])
    seed_rows: list[tuple[str, str, str, int, int]] = []
    for pointer_id, locus in sorted(registry.loci.items()):
        left = list(iter_syncmers(locus.left_flank, k=k, s=s, t=t))
        right = list(iter_syncmers(locus.right_flank, k=k, s=s, t=t))
        path_length = len(locus.left_flank) + 1 + len(locus.right_flank)
        conn.execute(
            "INSERT INTO loci VALUES(?,?,?,?,?,?,?,?,?,?)",
            (
                pointer_id,
                locus.chrom,
                locus.start,
                locus.end,
                locus.motif,
                locus.anchor_length,
                path_length,
                len(left) + len(right),
                len(left),
                len(right),
            ),
        )
        for hit in left:
            seed_rows.append((hit.seed, pointer_id, "L", hit.position, hit.orientation))
        right_offset = len(locus.left_flank) + 1
        for hit in right:
            seed_rows.append((hit.seed, pointer_id, "R", right_offset + hit.position, hit.orientation))
        if len(seed_rows) >= 100_000:
            conn.executemany("INSERT INTO seeds VALUES(?,?,?,?,?)", seed_rows)
            seed_rows.clear()
    if seed_rows:
        conn.executemany("INSERT INTO seeds VALUES(?,?,?,?,?)", seed_rows)
    conn.commit()
    conn.execute("ANALYZE")
    conn.close()
