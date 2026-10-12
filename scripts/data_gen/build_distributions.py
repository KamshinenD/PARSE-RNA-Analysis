"""Build per-(bp_type, edge_type, parameter) reference distributions.

Reads the three feature tables produced by extract_pair_features.py and emits

    data/reference/parameter_distributions.json

For each parameter consumed by the runtime scorer, computes per-cell:
  - weighted_median       M̂         weighted by 1 / (same-PDB count in cell)
  - weighted_d1, d9                  10th / 90th percentiles
  - scale_factor_upper    (d9 − M̂) / Φ⁻¹(0.9)        - Z' upper scale
  - scale_factor_lower    (M̂ − d1) / |Φ⁻¹(0.1)|      - Z' lower scale
  - n_observations, n_effective, n_pdbs

The output JSON includes fallback pools so the runtime can degrade to
broader data when a specific cell is sparse:
  bp_type/edge_type    primary cell        e.g. "G-C/cWW"
  _ANY/edge_type       pooled across bp_types
  _ANY/_ANY            global pool

Reference: Černý et al. 2026, NAR gkaf1335 - defines Z'.
"""

import argparse
import datetime
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _vendor.canonical import (  # noqa: E402
    SELF_RECIPROCAL_LW, SIGN_FLIPPING_PAIR_PARAMS,
)
DEFAULT_FEATURES_DIR = REPO_ROOT / "data" / "reference" / "high_quality_features"
DEFAULT_OUTPUT = (REPO_ROOT / "data" / "reference" / "scoring_tables"
                  / "parameter_distributions.json")

# Z' transform constants (paper Eq. 1)
PROBIT_90 = float(stats.norm.ppf(0.9))   # ≈ +1.2816
ABS_PROBIT_10 = abs(float(stats.norm.ppf(0.1)))

# Parameters that get Z' tables - must mirror what the runtime scorer reads.
PARAMETERS = {
    "pairs":  ["shear", "stretch", "stagger", "buckle", "propeller", "opening"],
    "hbonds": ["distance", "donor_angle", "acceptor_angle"],
}

# Categorical columns that get per-cell weighted frequency tables.
CATEGORICAL_COLUMNS = ["conformer"]


def canonicalize_dataframe(df: pd.DataFrame, source: str) -> pd.DataFrame:
    """Rewrite bp_type into canonical (alphabetical) order for self-reciprocal
    LW classes, flipping the sign of orientation-dependent pair-geometry
    parameters when a swap occurs.

    Returns a new DataFrame; the input is not mutated. ``source`` selects
    which sign-flip set applies - only "pairs" carries shear/stagger/buckle.
    """
    if "lw_class" not in df.columns or "bp_type" not in df.columns:
        return df
    df = df.copy()
    bp = df["bp_type"].astype(str)
    parts = bp.str.split("-", n=1, expand=True)
    # Rows that should be swapped: self-reciprocal LW + bp has a dash + left > right
    needs_swap = (
        df["lw_class"].isin(SELF_RECIPROCAL_LW)
        & bp.str.contains("-", regex=False)
        & (parts[0] > parts[1])
    )
    n = int(needs_swap.sum())
    if n == 0:
        return df
    print(f"  [canonicalize/{source}] swapped {n:,} rows into canonical orientation")
    df.loc[needs_swap, "bp_type"] = parts.loc[needs_swap, 1] + "-" + parts.loc[needs_swap, 0]
    if source == "pairs":
        for p in SIGN_FLIPPING_PAIR_PARAMS:
            if p in df.columns:
                df.loc[needs_swap, p] = -df.loc[needs_swap, p]
    return df


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    """Weighted quantile via cumulative-weight interpolation.

    Args:
        values: 1-D array of observations.
        weights: matching weights (positive).
        q: quantile in [0, 1].

    Returns:
        The interpolated value at quantile q. Returns nan for empty input.
    """
    if values.size == 0:
        return float("nan")
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cum = np.cumsum(w)
    total = cum[-1]
    if total <= 0:
        return float("nan")
    # Cumulative midpoints: paper-style "weighted percentile" uses (cum - w/2) / total
    p = (cum - w / 2.0) / total
    # Clamp q to the supported range, then linear-interpolate
    return float(np.interp(q, p, v))


def per_pdb_weights(pdb_ids: pd.Series) -> np.ndarray:
    """Per-PDB downweighting: each obs gets weight 1 / (count of same-PDB obs in this cell).

    Sum of weights across all obs from one PDB == 1, so a PDB contributes at
    most one effective sample regardless of how many pairs it has in the cell.
    Used by Černý et al. 2026 to fight ribosome / homolog dominance.
    """
    counts = pdb_ids.value_counts()
    return np.array([1.0 / counts[pdb_id] for pdb_id in pdb_ids])


def cell_distribution(values: np.ndarray, weights: np.ndarray, pdb_ids: pd.Series) -> dict:
    """Compute Z' parameters for one cell × one parameter."""
    mask = ~np.isnan(values)
    values = values[mask]
    weights = weights[mask]
    pdb_ids_masked = pdb_ids[mask]

    base = {
        "n_observations": int(values.size),
        "n_effective": round(float(weights.sum()), 2) if values.size else 0.0,
        "n_pdbs": int(pdb_ids_masked.nunique()),
        "weighted_median": None,
        "scale_factor_upper": None,
        "scale_factor_lower": None,
    }
    if values.size < 3:
        return base

    median = weighted_quantile(values, weights, 0.5)
    d9 = weighted_quantile(values, weights, 0.9)
    d1 = weighted_quantile(values, weights, 0.1)
    scale_upper = (d9 - median) / PROBIT_90
    scale_lower = (median - d1) / ABS_PROBIT_10

    base["weighted_median"] = round(median, 4)
    base["weighted_d1"] = round(d1, 4)
    base["weighted_d9"] = round(d9, 4)
    # Degenerate cells (zero spread): scales stay None so the scorer skips them.
    if scale_upper > 0 and scale_lower > 0:
        base["scale_factor_upper"] = round(scale_upper, 5)
        base["scale_factor_lower"] = round(scale_lower, 5)
    return base


def categorical_cell(values: pd.Series, weights: np.ndarray) -> dict:
    """Compute weighted frequency table for a categorical column in one cell.

    Returns
    -------
    dict with::
        n_observations:     raw count (excluding NaN)
        n_effective:        Σ weights
        n_pdbs:             distinct PDB count (computed by caller via .name path)
        frequencies:        {category_value: weighted_fraction}, sorted desc
    """
    mask = values.notna()
    v = values[mask].to_numpy()
    w = weights[mask]
    if v.size == 0:
        return {
            "n_observations": 0,
            "n_effective": 0.0,
            "frequencies": {},
        }
    total_weight = float(w.sum())
    freqs: dict[str, float] = {}
    for cat in pd.unique(v):
        cat_weight = float(w[v == cat].sum())
        freqs[str(cat)] = round(cat_weight / total_weight, 4)
    # Sort by frequency descending so the dominant category is first.
    freqs = dict(sorted(freqs.items(), key=lambda kv: -kv[1]))
    return {
        "n_observations": int(v.size),
        "n_effective": round(total_weight, 2),
        "frequencies": freqs,
    }


def build_categorical_for_dataframe(df: pd.DataFrame, columns: list[str]) -> dict:
    """Build per-(bp_type, lw_class) weighted frequency tables.

    Mirrors the cell + fallback layout used by build_for_dataframe so the
    runtime scorer walks the same hierarchy.
    """
    out: dict[str, dict[str, dict]] = {c: {} for c in columns}

    # Primary cells: every observed (bp_type, lw_class) combination
    for (bp_type, lw_class), group in df.groupby(["bp_type", "lw_class"]):
        weights = per_pdb_weights(group["pdb_id"])
        n_pdbs = group["pdb_id"].nunique()
        key = f"{bp_type}/{lw_class}"
        for c in columns:
            entry = categorical_cell(group[c], weights)
            entry["n_pdbs"] = int(n_pdbs)
            out[c][key] = entry

    # Edge-only pool
    for lw_class, group in df.groupby("lw_class"):
        weights = per_pdb_weights(group["pdb_id"])
        n_pdbs = group["pdb_id"].nunique()
        key = f"_ANY/{lw_class}"
        for c in columns:
            entry = categorical_cell(group[c], weights)
            entry["n_pdbs"] = int(n_pdbs)
            out[c][key] = entry

    # Global pool
    weights = per_pdb_weights(df["pdb_id"])
    n_pdbs = df["pdb_id"].nunique()
    for c in columns:
        entry = categorical_cell(df[c], weights)
        entry["n_pdbs"] = int(n_pdbs)
        out[c]["_ANY/_ANY"] = entry

    return out


def build_for_dataframe(df: pd.DataFrame, parameters: list[str]) -> dict:
    """Build the cell × parameter distribution table for one DataFrame.

    Returns a nested dict::
        {parameter_name: {"bp_type/edge_type": {z_prime_fields...}, ...}}

    Includes _ANY/edge_type pools and a global _ANY/_ANY pool for fallback.
    """
    out: dict[str, dict[str, dict]] = {p: {} for p in parameters}

    # Primary cells: every observed (bp_type, lw_class) combination
    for (bp_type, lw_class), group in df.groupby(["bp_type", "lw_class"]):
        weights = per_pdb_weights(group["pdb_id"])
        key = f"{bp_type}/{lw_class}"
        for p in parameters:
            out[p][key] = cell_distribution(
                group[p].to_numpy(dtype=float),
                weights,
                group["pdb_id"].reset_index(drop=True),
            )

    # Edge-only pool: pool across all bp_types per LW class
    for lw_class, group in df.groupby("lw_class"):
        weights = per_pdb_weights(group["pdb_id"])
        key = f"_ANY/{lw_class}"
        for p in parameters:
            out[p][key] = cell_distribution(
                group[p].to_numpy(dtype=float),
                weights,
                group["pdb_id"].reset_index(drop=True),
            )

    # Global pool: everything
    weights = per_pdb_weights(df["pdb_id"])
    for p in parameters:
        out[p]["_ANY/_ANY"] = cell_distribution(
            df[p].to_numpy(dtype=float),
            weights,
            df["pdb_id"].reset_index(drop=True),
        )

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--features-dir", type=Path, default=DEFAULT_FEATURES_DIR)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = ap.parse_args()

    print(f"Reading feature tables from {args.features_dir}/")
    pairs_df = pd.read_csv(args.features_dir / "pairs.csv")
    hbonds_df = pd.read_csv(args.features_dir / "hbonds.csv")
    torsions_df = pd.read_csv(args.features_dir / "torsions.csv")
    print(f"  pairs:    {len(pairs_df):,} rows")
    print(f"  hbonds:   {len(hbonds_df):,} rows")
    print(f"  torsions: {len(torsions_df):,} rows")

    print("\nCanonicalizing self-reciprocal LW orientations...")
    pairs_df = canonicalize_dataframe(pairs_df, "pairs")
    hbonds_df = canonicalize_dataframe(hbonds_df, "hbonds")
    torsions_df = canonicalize_dataframe(torsions_df, "torsions")

    print("\nBuilding Z' distributions for continuous parameters...")
    z_prime_distributions = {
        "pairs": build_for_dataframe(pairs_df, PARAMETERS["pairs"]),
        "hbonds": build_for_dataframe(hbonds_df, PARAMETERS["hbonds"]),
    }
    print("Building categorical frequency tables for backbone...")
    categorical_distributions = {
        "torsions": build_categorical_for_dataframe(torsions_df, CATEGORICAL_COLUMNS),
        # H-bond count is integer-valued and discrete: treat as categorical
        # so the runtime scorer can flag pairs whose count is in the rare
        # tail for their cell.
        "pairs": build_categorical_for_dataframe(pairs_df, ["num_base_hbonds"]),
    }

    output = {
        "generated_at": datetime.date.today().isoformat(),
        "method": {
            "z_prime": (
                "Z'(x) = (x − M̂) / scale_upper if x ≥ M̂, else / scale_lower. "
                "scale_upper = (D̂9 − M̂) / Φ⁻¹(0.9); scale_lower = (M̂ − D̂1) / |Φ⁻¹(0.1)|. "
                "Per Černý et al. 2026 (NAR gkaf1335)."
            ),
            "weighting": (
                "Per-PDB downweighting: each observation weight = 1 / (same-PDB "
                "count in cell). n_effective ≤ number of distinct PDBs in the cell."
            ),
            "fallback_hierarchy": [
                "<bp_type>/<edge_type>", "_ANY/<edge_type>", "_ANY/_ANY",
            ],
        },
        "z_prime_distributions": z_prime_distributions,
        "categorical_distributions": categorical_distributions,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    # Summary stats
    print(f"\nWrote {args.output}")
    for source, params in z_prime_distributions.items():
        n_cells = len(next(iter(params.values()))) if params else 0
        print(f"  Z'         {source:9s}: {len(params)} params × {n_cells} cells")
    for source, params in categorical_distributions.items():
        n_cells = len(next(iter(params.values()))) if params else 0
        print(f"  cat. freq. {source:9s}: {len(params)} columns × {n_cells} cells")

if __name__ == "__main__":
    main()


# Usage:
#   python3.12 scripts/reference/build_distributions.py
#   python3.12 scripts/reference/build_distributions.py \
#       --features-dir data/reference/features_strict \
#       --output       data/reference/parameter_distributions_strict.json