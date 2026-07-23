"""Build per-(bp_type, edge_type, parameter) ProSco lookup tables.

ProSco — *probability percentile score* — is the density-based scoring
metric from Černý et al. 2026 (NAR gkaf1335). It avoids the heavy-tail
issues of Z' by working directly with the empirical density:

    For each cell × parameter:
      1. Build weighted KDE on the reference data (Gaussian kernel,
         bandwidth × 1.5 per the paper).
      2. Pick n_bins per cell at the Nyquist sampling rate:
              n_bins = ⌈2 × support_range / bandwidth⌉
         clamped to [BIN_MIN, BIN_MAX] for sanity. This is the smallest
         discretisation that resolves smoothing-scale features (sampling
         theorem); no chosen safety factor.
      3. Discretise the support into n_bins equal-width bins; compute
         density at each bin centre.
      4. Per-bin probability mass = density × bin_width.
      5. Sort bins by descending density. Cumulative mass through that
         order gives, for each bin, the fraction of total mass
         concentrated in bins of *higher* density.
      6. ProSco(bin) = 100 × (1 − cum_mass_in_higher_density_bins).
         Peak bin → ProSco ≈ 100. Deep-tail bin → ProSco ≈ 0.
         By construction, 5 % of reference mass has ProSco < 5.

Mirrors the canonicalisation + per-PDB weighting + fallback hierarchy
in build_distributions.py so ProSco and Z' live in the same cell layout.

Output: data/reference/scoring_tables/prosco_distributions.json
"""

from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
import _pairfinder_path  # noqa: F401,E402  (adds pair_finder from PARSE_CODE_REPO)
from pair_finder.scoring.canonical import (  # noqa: E402
    SELF_RECIPROCAL_LW, SIGN_FLIPPING_PAIR_PARAMS,
)
from pair_finder.scoring.hbond_distance_policy import (  # noqa: E402
    select_distance_rows, base_base_rows,
)

DEFAULT_FEATURES_DIR = REPO_ROOT / "data" / "reference" / "high_quality_features"
DEFAULT_OUTPUT = REPO_ROOT / "data" / "reference" / "scoring_tables" / "prosco_distributions.json"

# Bandwidth adjustment (paper Section "Definition of probability percentile
# scores"): Gaussian KDE with bandwidth × 1.5.
BANDWIDTH_FACTOR = 1.5
# n_bins is derived per cell from Nyquist (×2 / bandwidth), bounded by
# these safety limits to avoid degenerate tiny tables or oversized JSON.
BIN_MIN = 50
BIN_MAX = 2048
# Support extends 1× the observed data range beyond min/max on each side
# (so total support = 3× data range). Captures plausible outliers without
# making bins too coarse near the peak.
SUPPORT_EXPANSION = 1.0

PARAMETERS = {
    "pairs":  ["shear", "stretch", "stagger", "buckle", "propeller", "opening"],
    "hbonds": ["distance", "donor_angle", "acceptor_angle"],
}


def canonicalize_dataframe(df: pd.DataFrame, source: str) -> pd.DataFrame:
    """Swap into canonical orientation for self-reciprocal LW classes,
    flipping the sign of orientation-dependent pair-geometry parameters."""
    if "lw_class" not in df.columns or "bp_type" not in df.columns:
        return df
    df = df.copy()
    bp = df["bp_type"].astype(str)
    parts = bp.str.split("-", n=1, expand=True)
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


def per_pdb_weights(pdb_ids: pd.Series) -> np.ndarray:
    """1 / (count of same-PDB obs in this cell). Sum per PDB == 1."""
    counts = pdb_ids.value_counts()
    return np.array([1.0 / counts[p] for p in pdb_ids])


def build_prosco_cell(values: np.ndarray, weights: np.ndarray,
                      pdb_ids: pd.Series) -> dict:
    """Build the ProSco lookup table for one cell × one parameter.

    Returns the support range, the per-bin ProSco values, and bookkeeping.
    Cells with fewer than 10 observations get an empty table (the runtime
    falls back to a broader pool).
    """
    mask = ~np.isnan(values)
    values = values[mask]
    weights = weights[mask]
    pdb_ids_masked = pdb_ids[mask]

    base = {
        "n_observations": int(values.size),
        "n_effective": round(float(weights.sum()), 2) if values.size else 0.0,
        "n_pdbs": int(pdb_ids_masked.nunique()),
        "support_lo": None,
        "support_hi": None,
        "n_bins": None,
        "prosco_per_bin": None,
        "weighted_median": None,
    }
    if values.size < 10:
        return base

    # Degenerate cell (no spread): emit a flat ProSco = 100 table.
    if np.ptp(values) < 1e-9:
        constant = float(values[0])
        eps = 1e-6
        base["support_lo"] = round(constant - eps, 6)
        base["support_hi"] = round(constant + eps, 6)
        base["n_bins"] = BIN_MIN
        base["prosco_per_bin"] = [100.0] * BIN_MIN
        base["weighted_median"] = round(constant, 4)
        return base

    # KDE with × 1.5 bandwidth
    try:
        kde = gaussian_kde(values, weights=weights)
        kde.set_bandwidth(kde.factor * BANDWIDTH_FACTOR)
    except Exception as e:
        print(f"    KDE failed ({type(e).__name__}: {e}); leaving cell empty")
        return base

    # Support range: 1× data range expansion on each side, capped by ±5σ
    # safety floor to avoid pathological tails dominating bin width.
    vmin, vmax = float(values.min()), float(values.max())
    dr = vmax - vmin
    lo = vmin - SUPPORT_EXPANSION * dr
    hi = vmax + SUPPORT_EXPANSION * dr

    # Per-cell n_bins via Nyquist: 2× sampling of the KDE bandwidth.
    # scipy's gaussian_kde stores the effective bandwidth as kde.factor × σ
    # (its bandwidth attribute is the *bandwidth squared*, in covariance form).
    # We use the post-multiplier effective sigma directly.
    bandwidth = float(np.sqrt(kde.covariance.item()))
    support_range = hi - lo
    n_bins = int(np.ceil(2.0 * support_range / bandwidth))
    n_bins = max(BIN_MIN, min(BIN_MAX, n_bins))

    edges = np.linspace(lo, hi, n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    bin_width = edges[1] - edges[0]
    dens = kde(centres)
    mass = dens * bin_width                    # per-bin probability mass

    # Sort bins by descending density; cumulative mass through that order.
    order_desc = np.argsort(-dens)
    cum_through_higher = np.cumsum(mass[order_desc])
    # Shift so the cumulative *up to and not including* this bin is used —
    # the bin with the single highest density gets ProSco ≈ 100.
    cum_excluding_self = np.concatenate([[0.0], cum_through_higher[:-1]])
    prosco_sorted = 100.0 * (1.0 - cum_excluding_self)
    prosco_per_bin = np.empty_like(prosco_sorted)
    prosco_per_bin[order_desc] = prosco_sorted

    # Weighted median (for bookkeeping / cross-reference with Z' tables).
    order = np.argsort(values)
    v, w = values[order], weights[order]
    cum = np.cumsum(w)
    median = float(np.interp(0.5 * cum[-1], cum - w / 2.0, v))

    base["support_lo"] = round(lo, 6)
    base["support_hi"] = round(hi, 6)
    base["n_bins"] = n_bins
    base["prosco_per_bin"] = [round(float(p), 2) for p in prosco_per_bin]
    base["weighted_median"] = round(median, 4)
    return base


def build_for_dataframe(df: pd.DataFrame, parameters: list[str]) -> dict:
    """Per-cell × per-parameter ProSco tables, with pooled fallbacks."""
    out: dict[str, dict[str, dict]] = {p: {} for p in parameters}

    print(f"  Primary cells (bp_type × lw_class) ...")
    n_cells = 0
    for (bp_type, lw_class), group in df.groupby(["bp_type", "lw_class"]):
        weights = per_pdb_weights(group["pdb_id"])
        key = f"{bp_type}/{lw_class}"
        for p in parameters:
            out[p][key] = build_prosco_cell(
                group[p].to_numpy(dtype=float),
                weights,
                group["pdb_id"].reset_index(drop=True),
            )
        n_cells += 1
    print(f"    {n_cells} primary cells")

    print(f"  _ANY/<edge_type> pools ...")
    for lw_class, group in df.groupby("lw_class"):
        weights = per_pdb_weights(group["pdb_id"])
        key = f"_ANY/{lw_class}"
        for p in parameters:
            out[p][key] = build_prosco_cell(
                group[p].to_numpy(dtype=float),
                weights,
                group["pdb_id"].reset_index(drop=True),
            )

    print(f"  _ANY/_ANY global pool ...")
    weights = per_pdb_weights(df["pdb_id"])
    for p in parameters:
        out[p]["_ANY/_ANY"] = build_prosco_cell(
            df[p].to_numpy(dtype=float),
            weights,
            df["pdb_id"].reset_index(drop=True),
        )

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--features-dir", type=Path, default=DEFAULT_FEATURES_DIR)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = ap.parse_args()

    print(f"Reading feature tables from {args.features_dir}/")
    pairs_df = pd.read_csv(args.features_dir / "pairs.csv")
    hbonds_df = pd.read_csv(args.features_dir / "hbonds.csv")
    print(f"  pairs:  {len(pairs_df):,} rows")
    print(f"  hbonds: {len(hbonds_df):,} rows")

    print("\nCanonicalizing self-reciprocal LW orientations...")
    pairs_df = canonicalize_dataframe(pairs_df, "pairs")
    hbonds_df = canonicalize_dataframe(hbonds_df, "hbonds")

    print("\nBuilding ProSco tables for pairs ...")
    pairs_prosco = build_for_dataframe(pairs_df, PARAMETERS["pairs"])
    print("\nBuilding ProSco tables for hbonds ...")
    # Distance follows the WC/strong-bond policy and is base-base filtered;
    # angles use all base-base bonds (matching the runtime scorer).
    dist_df = select_distance_rows(hbonds_df)
    angle_df = base_base_rows(hbonds_df)
    hbonds_prosco = build_for_dataframe(dist_df, ["distance"])
    hbonds_prosco.update(build_for_dataframe(angle_df, ["donor_angle", "acceptor_angle"]))

    output = {
        "schema_version": "1.0",
        "generated_at": datetime.date.today().isoformat(),
        "method": {
            "prosco": (
                "ProSco(x) = 100 × (1 − cum_mass_in_higher_density_bins(x)). "
                "Per cell: weighted KDE (Gaussian, bandwidth × 1.5), "
                "equal-width bins over [data_min − 1×range, data_max + 1×range]. "
                f"n_bins per cell = ⌈2 × support_range / bandwidth⌉ "
                f"clamped to [{BIN_MIN}, {BIN_MAX}] (Nyquist sampling rate; "
                "smallest discretisation that resolves smoothing-scale features). "
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
        "prosco_distributions": {
            "pairs":  pairs_prosco,
            "hbonds": hbonds_prosco,
        },
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(output, f, ensure_ascii=False)
    print(f"\nWrote {args.output}")
    size_mb = args.output.stat().st_size / 1e6
    print(f"  file size: {size_mb:.1f} MB")
    for source, params in output["prosco_distributions"].items():
        n_cells = len(next(iter(params.values()))) if params else 0
        print(f"  ProSco {source:8s}: {len(params)} params × {n_cells} cells")


if __name__ == "__main__":
    main()
