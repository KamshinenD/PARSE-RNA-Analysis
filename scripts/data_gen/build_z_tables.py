"""Build weighted Z' tables (median + 1st/9th deciles) from the GOOD reference.

Generates z_tables.json with per-parameter Z' statistics (M, sl, su) used by
the Cerny severity function for continuous parameter scoring.

Structure mirrors prosco_distributions.json:
{
  "pairs": { "shear": { "G-C/cWW": {M, sl, su, n}, "_ANY/cWW": ..., "_ANY/_ANY": ... }, ... },
  "hbonds": { "distance": { ... }, ... }
}

M = weighted median
sl = sigma_lower = (M - D1) / 1.2816  (D1 = 10th percentile)
su = sigma_upper = (D9 - M) / 1.2816  (D9 = 90th percentile)

This is the first step in building a new RNA quality scorer. Run this to
generate fresh Z' tables from data/reference/high_quality_features/ before running
calculate_penalty_weights.py.
"""
import json
import sys
import numpy as np
import pandas as pd
from pathlib import Path

# Repository and paths
REPO = Path(__file__).resolve().parent.parent.parent
import _pairfinder_path  # noqa: F401,E402  (adds pair_finder from PARSE_CODE_REPO)
from pair_finder.scoring.hbond_distance_policy import (  # noqa: E402
    select_distance_rows, base_base_rows,
)
REF = REPO / "data" / "reference"
OUT = REF / "scoring_tables" / "z_tables.json"

PHI_INV_09 = 1.2815515655446004  # Phi^-1(0.9) — same scaling as Cerny

PAIR_PARAMS = ["shear", "stretch", "stagger", "buckle", "propeller", "opening"]
HB_PARAMS = ["distance", "donor_angle", "acceptor_angle"]
SELF_RECIP = {"cWW", "tWW", "cHH", "tHH", "cSS", "tSS"}
SIGN_FLIP = {"shear", "stagger", "buckle"}


def per_pdb_weights(pdb_ids: pd.Series) -> np.ndarray:
    """Per-PDB weighting: each structure contributes equally."""
    counts = pdb_ids.value_counts()
    return np.array([1.0 / counts[p] for p in pdb_ids])


def weighted_quantile(values, weights, q):
    """Compute weighted quantile q ∈ [0, 1]."""
    order = np.argsort(values)
    v, w = values[order], weights[order]
    cum = np.cumsum(w)
    return float(np.interp(q * cum[-1], cum - w / 2.0, v))


def build_cell(values: np.ndarray, weights: np.ndarray):
    """Build Z' cell statistics from weighted values."""
    mask = ~np.isnan(values)
    v, w = values[mask], weights[mask]
    if len(v) < 10 or np.ptp(v) < 1e-9:
        return None
    M = weighted_quantile(v, w, 0.5)
    D1 = weighted_quantile(v, w, 0.1)
    D9 = weighted_quantile(v, w, 0.9)
    sl = (M - D1) / PHI_INV_09
    su = (D9 - M) / PHI_INV_09
    return {
        "M": round(M, 6),
        "sl": round(max(sl, 1e-9), 6),
        "su": round(max(su, 1e-9), 6),
        "n": int(len(v)),
    }


def canon_bp(bp, lw):
    """Canonicalize base-pair type for self-reciprocal LW classes."""
    if lw not in SELF_RECIP or "-" not in bp:
        return bp, False
    a, b = bp.split("-", 1)
    return (bp, False) if a <= b else (b + "-" + a, True)


def build_tables(df: pd.DataFrame, params: list[str], source: str) -> dict:
    """Build Z' tables for all parameters in df."""
    out = {p: {} for p in params}
    # Canonicalize + sign-flip for pair params
    d = df.copy()
    if source == "pairs":
        swap = d["lw_class"].isin(SELF_RECIP) & (
            d["bp_type"].str.split("-").str[0] > d["bp_type"].str.split("-").str[1]
        )
        for p in SIGN_FLIP:
            if p in d.columns:
                d.loc[swap, p] = -d.loc[swap, p]
    d["_cb"] = [canon_bp(b, l)[0] for b, l in zip(d["bp_type"], d["lw_class"])]

    def add(key, grp):
        w = per_pdb_weights(grp["pdb_id"])
        for p in params:
            if p not in grp.columns:
                continue
            c = build_cell(grp[p].to_numpy(float), w)
            if c:
                out[p][key] = c

    for (cb, lw), g in d.groupby(["_cb", "lw_class"]):
        add(f"{cb}/{lw}", g)
    for lw, g in d.groupby("lw_class"):
        add(f"_ANY/{lw}", g)
    add("_ANY/_ANY", d)
    return out


def main():
    print("Loading GOOD reference tables...")
    good_pairs = pd.read_csv(REF / "high_quality_features/pairs.csv")
    good_hb = pd.read_csv(REF / "high_quality_features/hbonds.csv")
    print(f"  pairs: {len(good_pairs)} rows / {good_pairs.pdb_id.nunique()} structs")
    print(f"  hbonds: {len(good_hb)} rows")

    print("\nBuilding Z' tables for pair parameters...")
    pair_tables = build_tables(good_pairs, PAIR_PARAMS, "pairs")
    print("Building Z' tables for H-bond parameters...")
    # Distance follows the WC/strong-bond policy (base-base filtered); angles
    # use all base-base bonds, matching the runtime scorer.
    hb_tables = build_tables(select_distance_rows(good_hb), ["distance"], "hbonds")
    hb_tables.update(
        build_tables(base_base_rows(good_hb), ["donor_angle", "acceptor_angle"], "hbonds")
    )

    result = {"pairs": pair_tables, "hbonds": hb_tables}
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nWrote {OUT}")
    for src, tab in result.items():
        for p, cells in tab.items():
            n = sum(1 for v in cells.values() if v)
            print(f"  {src}/{p}: {n} cells")


if __name__ == "__main__":
    main()
