#!/usr/bin/env bash
set -euo pipefail
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_root"
python -m compileall -q src
python -m pytest -q
python -m strpg.cli --help >/dev/null
python -m strpg.cli genotype --help >/dev/null
bash scripts/run_demo.sh
bash scripts/run_phmm_demo.sh
