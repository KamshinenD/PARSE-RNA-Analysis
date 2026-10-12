"""Run the PARSE C++ engine (`parse --details`) and return its JSON.

The engine is the single source of features: this replaces the old Python
pair-finder pipeline. Build the `parse` binary from the PARSE-RNA repo
(https://github.com/KamshinenD/PARSE-RNA) and either put it on PATH or point
PARSE_ENGINE at it. Structures are read from data/structures/<ID>.cif
(scripts/fetch_structures.py); an entry that is not there is downloaded by the
engine into data/structures/engine_cache/, so nothing is read or written outside
the repository.
"""

import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
STRUCTURES = REPO / "data" / "structures"
ENGINE_CACHE = STRUCTURES / "engine_cache"


def structure_args(pdb_id: str) -> list[str]:
    """Engine input for one entry: the local mmCIF if present, else the id with the
    download cache inside data/structures."""
    local = STRUCTURES / f"{pdb_id.upper()}.cif"
    if local.exists():
        return [str(local), "--no-download"]
    return [pdb_id, "--cache-dir", str(ENGINE_CACHE)]


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
        src = [pdb_id, "--cache-dir", str(cache_dir)] if cache_dir else structure_args(pdb_id)
        cmd = [engine, *src, "--details", "--out", out]
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


# DNA residue codes. Mirrors BaseTyping::is_dna in the engine.
_DNA_LITERAL = frozenset({"DA", "DC", "DG", "DT", "DU"})


def is_dna_name(name: str) -> bool:
    """True when a residue name is DNA."""
    n = (name or "").strip().upper()
    return n in _DNA_LITERAL


def is_dna_res_id(res_id: str) -> bool:
    """True when the residue in a '<chain>-<name>-<seq>' id is DNA."""
    return is_dna_name(res_name(res_id))


def is_scorable(pair: dict) -> bool:
    """base-base, not ambiguous, single (no '|') lw_class, neither side DNA.

    Tables are calibrated on RNA, so either side DNA drops the pair, hybrids
    included. Must match Scorer::is_scorable in the engine.
    """
    c = pair.get("classification") or {}
    lw = c.get("lw_class") or ""
    if pair.get("pair_category") != "base-base" or c.get("is_ambiguous"):
        return False
    if not lw or "|" in lw:
        return False
    if is_dna_res_id(pair.get("res_id1", "")) or is_dna_res_id(pair.get("res_id2", "")):
        return False
    return True


def res_name(res_id: str) -> str:
    """Raw residue name (middle field of '<chain>-<name>-<seq>')."""
    return res_id.split("-")[1]
