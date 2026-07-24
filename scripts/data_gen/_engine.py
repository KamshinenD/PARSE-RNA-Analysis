"""Run the PARSE C++ engine (`parse --details`) and return its JSON.

The engine is the single source of features: this replaces the old Python
pair-finder pipeline. Build the `parse` binary from the PARSE-RNA repo
(https://github.com/KamshinenD/PARSE-RNA) and either put it on PATH or point
PARSE_ENGINE at it. `parse` downloads the mmCIF for each PDB id itself (cached
under ~/.cache/parse), so no local structure files are needed.
"""

import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


def find_engine(explicit: str | None = None) -> str:
    for cand in (explicit, os.environ.get("PARSE_ENGINE"), shutil.which("parse")):
        if cand and Path(cand).exists():
            return str(Path(cand).resolve())
    raise FileNotFoundError(
        "PARSE engine not found. Build PARSE-RNA and set PARSE_ENGINE=/path/to/parse "
        "(or put `parse` on PATH)."
    )


def run_parse(pdb_id: str, engine: str, cache_dir: str | None = None) -> dict:
    """Return the `parse <pdb_id> --details` JSON (downloads the mmCIF if needed)."""
    fd, out = tempfile.mkstemp(suffix=".json")
    os.close(fd)
    try:
        cmd = [engine, pdb_id, "--details", "--out", out]
        if cache_dir:
            cmd += ["--cache-dir", str(cache_dir)]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0 or os.path.getsize(out) == 0:
            raise RuntimeError(f"parse failed for {pdb_id}: {r.stderr.strip()[:200]}")
        with open(out) as f:
            return json.load(f)
    finally:
        try:
            os.unlink(out)
        except OSError:
            pass


def is_scorable(pair: dict) -> bool:
    """base-base, not ambiguous, single (no '|') lw_class — the scored pairs."""
    c = pair.get("classification") or {}
    lw = c.get("lw_class") or ""
    return pair.get("pair_category") == "base-base" and not c.get("is_ambiguous") and lw and "|" not in lw


def res_name(res_id: str) -> str:
    """Raw residue name (middle field of '<chain>-<name>-<seq>')."""
    return res_id.split("-")[1]
