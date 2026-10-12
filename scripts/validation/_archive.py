"""Loaders for the engine-output archive the re-solved cohort is built from.

Everything is read from data/: the per-residue and per-pair engine output over the
non-redundant set and over the extra members of each BGSU equivalence class
(data/reference/engine_archive/), and the PDB metadata.

    residues_orig     one row per classifiable residue of every non-redundant entry:
                      pdb_id, res_id, chain, num, base_type, suiteness, is_outlier,
                      pucker_outlier, pperp, s_delta
    residues_expand   the same (fewer columns) for the extra class members
    pairs_orig        one row per selected base pair: pdb_id, res_id1, res_id2,
                      lw_class, score, tier
    pairs_expand      the same for the extra class members
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
ARCHIVE = DATA / "reference" / "engine_archive"
STRUCTURES = DATA / "structures"


def residues_orig() -> pd.DataFrame:
    return pd.read_parquet(ARCHIVE / "residues_orig.parquet")


def residues_expand() -> pd.DataFrame:
    return pd.read_parquet(ARCHIVE / "residues_expand.parquet")


def pairs_orig() -> pd.DataFrame:
    return pd.read_parquet(ARCHIVE / "pairs_orig.parquet")


def pairs_expand() -> pd.DataFrame:
    return pd.read_parquet(ARCHIVE / "pairs_expand.parquet")


def uniques() -> pd.DataFrame:
    return pd.read_csv(DATA / "reference" / "reference_sets" / "uniqueRNAS.csv")


def metadata() -> pd.DataFrame:
    """pdb_id, resolution, method and nucleotide count from metadata_cache.json."""
    md = json.loads((DATA / "reference" / "metadata" / "metadata_cache.json").read_text())
    nt = md.get("nucleotide_counts", {})
    rows = []
    for pdb, m in md["validation_metrics"].items():
        if not isinstance(m, dict):
            continue
        res = m.get("Refinement_resolution")
        if res is None:
            res = m.get("EM Resolution (Å)")
        rows.append({"pdb_id": pdb, "resolution": res,
                     "method": m.get("Experimental_method"), "n_nt": nt.get(pdb)})
    df = pd.DataFrame(rows)
    df["resolution"] = pd.to_numeric(df.resolution, errors="coerce")
    df["n_nt"] = pd.to_numeric(df.n_nt, errors="coerce")
    return df
