"""Put the PARSE reference code (pair_finder) on sys.path for the generation scripts.

These scripts regenerate the scoring tables (ProSco / Z' / penalty weights) from the
shipped feature tables, so they need the scoring code. The figure notebooks do NOT
need this — only scripts/data_gen/ does. Point PARSE_CODE_REPO at the pair_finder
repository root (the directory containing ``src/pair_finder``):

    export PARSE_CODE_REPO=/path/to/pair-finder
"""
import os
import sys
from pathlib import Path

_code = os.environ.get("PARSE_CODE_REPO")
if not _code:
    raise SystemExit(
        "PARSE_CODE_REPO is not set. Point it at the pair_finder repository root "
        "(the directory containing src/pair_finder) to run the scoring-table "
        "generation scripts. The figure notebooks do not need this.")
_src = Path(_code).expanduser() / "src"
if not (_src / "pair_finder").is_dir():
    raise SystemExit(f"PARSE_CODE_REPO={_code!r} has no src/pair_finder.")
sys.path.insert(0, str(_src))
