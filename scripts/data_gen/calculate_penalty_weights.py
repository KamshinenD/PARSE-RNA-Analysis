"""Derive per-parameter penalty weights using Shannon (log-odds) approach.

For each of the 9 parameters individually:
  - Compute mean Cerny severity in GOOD set vs POOR set
  - weight_i = log(s_poor_i / s_good_i)  [clipped to >= 0]
  - Renormalize to sum = 100

This is the pure information-theoretic approach. Per-parameter weights are
derived independently, then renormalized to sum to 100 for the penalty
calculation.

Output: data/reference/scoring_tables/penalty_weights.json

This script must stay in sync with src/pair_finder/scoring/issues.py. The
severity calculation here mirrors the Cerny function used at runtime.

Reference: Cerny et al. NAR 2026 (gkaf1335). Severity = 0 if ProSco >= 5
(Preferred), else min(1, |Z'|/5).
"""

import argparse
import datetime
import json
import math
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_FEATURES_DIR = REPO_ROOT / "data" / "reference" / "high_quality_features"
DEFAULT_POOR_FEATURES_DIR = REPO_ROOT / "data" / "reference" / "low_quality_features"
DEFAULT_DISTS = REPO_ROOT / "data" / "reference" / "scoring_tables" / "parameter_distributions.json"
DEFAULT_PROSCO_DISTS = REPO_ROOT / "data" / "reference" / "scoring_tables" / "prosco_distributions.json"
DEFAULT_Z_TABLES = REPO_ROOT / "data" / "reference" / "scoring_tables" / "z_tables.json"
DEFAULT_OUTPUT = REPO_ROOT / "data" / "reference" / "scoring_tables" / "penalty_weights.json"

from _vendor.canonical import canonicalize_bp, maybe_flip  # noqa: E402
from _vendor.hbond_distance_policy import distances_to_score  # noqa: E402
from _vendor.helpers import REPORT_GROUPS  # noqa: E402

# ---------------------------------------------------------------------------
# Parameters and groupings
# ---------------------------------------------------------------------------

PAIR_PARAMS = ["shear", "stretch", "stagger", "buckle", "propeller", "opening"]
HB_PARAMS = ["distance", "hbond_angles"]  # Combined donor + acceptor angles
ALL_PARAMS = PAIR_PARAMS + HB_PARAMS  # 8 continuous parameters
FEATURES = ALL_PARAMS + ["incorrect_hbond_count"]  # 9 total

SELF_RECIP = {"cWW", "tWW", "cHH", "tHH", "cSS", "tSS"}
SIGN_FLIP = {"shear", "stagger", "buckle"}
BASE_ATOMS = {"N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"}

# ---------------------------------------------------------------------------
# ProSco and Z' lookup
# ---------------------------------------------------------------------------


class ProScoLookup:
    """ProSco lookup with bp/lw → _ANY/lw → _ANY/_ANY fallback."""

    def __init__(self, prosco_distributions: dict):
        self.dist = prosco_distributions.get(
            "prosco_distributions", prosco_distributions
        )

    def _cell(self, source: str, param: str, bp_type: str, lw_class: str):
        canonical_bp, _ = canonicalize_bp(bp_type, lw_class)
        table = self.dist.get(source, {}).get(param, {})
        for key in (f"{canonical_bp}/{lw_class}", f"_ANY/{lw_class}", "_ANY/_ANY"):
            entry = table.get(key)
            if entry and entry.get("prosco_per_bin"):
                return entry
        return None

    def prosco(self, value, source: str, param: str, bp_type: str, lw_class: str):
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return None
        _, swapped = canonicalize_bp(bp_type, lw_class)
        value = maybe_flip(value, source, param, swapped)
        cell = self._cell(source, param, bp_type, lw_class)
        if cell is None:
            return None
        lo, hi = cell["support_lo"], cell["support_hi"]
        n_bins = cell["n_bins"]
        idx = int((value - lo) / (hi - lo) * n_bins)
        idx = max(0, min(n_bins - 1, idx))
        return float(cell["prosco_per_bin"][idx])


class ZPrimeLookup:
    """Z' lookup with bp/lw → _ANY/lw → _ANY/_ANY fallback."""

    def __init__(self, z_tables: dict):
        self.tables = z_tables

    def _cell(self, source: str, param: str, bp_type: str, lw_class: str):
        canonical_bp, _ = canonicalize_bp(bp_type, lw_class)
        table = self.tables.get(source, {}).get(param, {})
        for key in (f"{canonical_bp}/{lw_class}", f"_ANY/{lw_class}", "_ANY/_ANY"):
            entry = table.get(key)
            if entry:
                return entry
        return None

    def zprime(self, value, source: str, param: str, bp_type: str, lw_class: str):
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return None
        _, swapped = canonicalize_bp(bp_type, lw_class)
        value = maybe_flip(value, source, param, swapped)
        cell = self._cell(source, param, bp_type, lw_class)
        if cell is None:
            return None
        M, sl, su = cell["M"], cell["sl"], cell["su"]
        if value >= M:
            return (value - M) / su if su > 1e-9 else None
        else:
            return (value - M) / sl if sl > 1e-9 else None


# ---------------------------------------------------------------------------
# Cerny severity calculation
# ---------------------------------------------------------------------------


def cerny_severity(
    value, prosco_lookup, zprime_lookup, source, param, bp_type, lw_class
):
    """Cerny severity: 0 if ProSco >= 5 (Preferred), else min(1, |Z'|/5)."""
    p = prosco_lookup.prosco(value, source, param, bp_type, lw_class)
    if p is None:
        return 0.0
    if p >= 5.0:
        return 0.0
    z = zprime_lookup.zprime(value, source, param, bp_type, lw_class)
    if z is None:
        return 1.0
    return min(1.0, abs(z) / 5.0)


def canonical_count(cat_dists, bp_type, lw_class):
    """Modal num_base_hbonds for (bp_type, lw_class) cell."""
    canonical_bp, _ = canonicalize_bp(bp_type, lw_class)
    table = cat_dists.get("pairs", {}).get("num_base_hbonds", {})
    for key in (f"{canonical_bp}/{lw_class}", f"_ANY/{lw_class}", "_ANY/_ANY"):
        entry = table.get(key, {})
        f = entry.get("frequencies") if isinstance(entry, dict) else None
        if f:
            return int(max(f, key=lambda kk: f[kk]))
    return 1


# ---------------------------------------------------------------------------
# Compute severities for a pair dataset
# ---------------------------------------------------------------------------


def compute_severities(
    pairs_df, hbonds_df, prosco_lookup, zprime_lookup, cat_dists
):
    """Returns DataFrame with one severity column per parameter + pdb_id."""
    hb = hbonds_df[
        hbonds_df.donor_atom.isin(BASE_ATOMS)
        & hbonds_df.acceptor_atom.isin(BASE_ATOMS)
    ]
    hbg = {k: v for k, v in hb.groupby(["pdb_id", "pair_idx"])}

    rows = []
    for r in pairs_df.itertuples(index=False):
        bp, lw = r.bp_type, r.lw_class
        _, sw = canonicalize_bp(bp, lw)

        def pv(p):
            v = getattr(r, p)
            return -v if (sw and p in SIGN_FLIP) else v

        row = {"pdb_id": r.pdb_id}
        # 6 pair geometry parameters individually
        for p in PAIR_PARAMS:
            row[p] = cerny_severity(pv(p), prosco_lookup, zprime_lookup, "pairs", p, bp, lw)

        # H-bond parameters across base-base bonds. Distance follows the
        # WC/strong-bond policy (mirrors issues.py); angles take max over bonds.
        sub = hbg.get((r.pdb_id, r.pair_idx))
        d_sev = angle_sev = 0.0
        if sub is not None:
            bonds = list(sub.itertuples(index=False))
            dist_targets = distances_to_score([h.distance for h in bonds], bp, lw)
            d_sev = max(
                (cerny_severity(x, prosco_lookup, zprime_lookup, "hbonds", "distance", bp, lw)
                 for x in dist_targets),
                default=0.0,
            )
            for h in bonds:
                # Combine donor and acceptor angles into single hbond_angles severity
                donor_sev = cerny_severity(
                    h.donor_angle, prosco_lookup, zprime_lookup, "hbonds", "donor_angle", bp, lw
                )
                accept_sev = cerny_severity(
                    h.acceptor_angle, prosco_lookup, zprime_lookup, "hbonds", "acceptor_angle", bp, lw
                )
                angle_sev = max(angle_sev, donor_sev, accept_sev)
        row["distance"] = d_sev
        row["hbond_angles"] = angle_sev

        # incorrect hbond count: (canonical - actual) / canonical
        can = canonical_count(cat_dists, bp, lw)
        act = int(r.num_base_hbonds)
        row["incorrect_hbond_count"] = (
            0.0 if (can <= 0 or act >= can) else (can - act) / can
        )

        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Derive per-parameter penalty weights using Shannon log-odds."
    )
    parser.add_argument(
        "--features-dir",
        type=Path,
        default=DEFAULT_FEATURES_DIR,
        help="Path to GOOD reference features/ directory",
    )
    parser.add_argument(
        "--poor-features-dir",
        type=Path,
        default=DEFAULT_POOR_FEATURES_DIR,
        help="Path to POOR reference low_quality_features/ directory",
    )
    parser.add_argument(
        "--distributions",
        type=Path,
        default=DEFAULT_DISTS,
        help="Path to parameter_distributions.json",
    )
    parser.add_argument(
        "--prosco-distributions",
        type=Path,
        default=DEFAULT_PROSCO_DISTS,
        help="Path to prosco_distributions.json",
    )
    parser.add_argument(
        "--z-tables",
        type=Path,
        default=DEFAULT_Z_TABLES,
        help="Path to z_tables.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Path to output penalty_weights.json",
    )
    args = parser.parse_args()

    # Load distributions
    print("Loading distributions...")
    with open(args.distributions) as f:
        param_dists = json.load(f)
    cat_dists = param_dists["categorical_distributions"]

    with open(args.prosco_distributions) as f:
        prosco_payload = json.load(f)
    prosco_lookup = ProScoLookup(prosco_payload)

    with open(args.z_tables) as f:
        z_tables = json.load(f)
    zprime_lookup = ZPrimeLookup(z_tables)

    # Load reference tables
    print("Loading reference tables...")
    good_p = pd.read_csv(args.features_dir / "pairs.csv")
    good_h = pd.read_csv(args.features_dir / "hbonds.csv")
    poor_p = pd.read_csv(args.poor_features_dir / "pairs.csv")
    poor_h = pd.read_csv(args.poor_features_dir / "hbonds.csv")
    print(f"  GOOD:  {len(good_p)} pairs / {good_p.pdb_id.nunique()} structs")
    print(f"  POOR:  {len(poor_p)} pairs / {poor_p.pdb_id.nunique()} structs")

    # Compute severities
    print("\nComputing severities for GOOD set...")
    g_sev = compute_severities(good_p, good_h, prosco_lookup, zprime_lookup, cat_dists)
    print("Computing severities for POOR set...")
    p_sev = compute_severities(poor_p, poor_h, prosco_lookup, zprime_lookup, cat_dists)

    # Compute mean severity for each parameter in each set
    print("\n" + "=" * 70)
    print("SHANNON (LOG-ODDS) WEIGHTS - Individual Parameters")
    print("=" * 70)

    n_good = len(g_sev)
    n_poor = len(p_sev)

    # Laplace smoothing: (sum + 1) / (n + 2) to avoid log(0)
    results = []
    for param in FEATURES:
        mean_good = g_sev[param].mean()
        mean_poor = p_sev[param].mean()

        # Laplace smoothing
        s_good = (mean_good * n_good + 1) / (n_good + 2)
        s_poor = (mean_poor * n_poor + 1) / (n_poor + 2)

        # Log-odds (clip negative to 0)
        if s_good > 0:
            raw_logodds = max(0.0, math.log(s_poor / s_good))
        else:
            raw_logodds = 0.0

        ratio = s_poor / s_good if s_good > 0 else 0.0

        results.append(
            {
                "param": param,
                "mean_good": mean_good,
                "mean_poor": mean_poor,
                "s_good": s_good,
                "s_poor": s_poor,
                "ratio": ratio,
                "raw_logodds": raw_logodds,
            }
        )

    # Renormalize to sum = 100
    total_logodds = sum(r["raw_logodds"] for r in results)
    if total_logodds < 1e-9:
        raise ValueError("All log-odds are zero or negative — something is wrong")

    for r in results:
        r["weight"] = 100.0 * r["raw_logodds"] / total_logodds

    # Sort by weight descending
    results.sort(key=lambda r: -r["weight"])

    print(
        f"\n{'Parameter':<25} {'Mean_GOOD':<10} {'Mean_POOR':<10} {'Ratio':<8} {'LogOdds':<10} {'Weight':<8}"
    )
    print("-" * 80)
    for r in results:
        print(
            f"{r['param']:<25} {r['mean_good']:<10.4f} {r['mean_poor']:<10.4f} "
            f"{r['ratio']:<8.2f} {r['raw_logodds']:<10.4f} {r['weight']:<8.2f}"
        )

    total_weight = sum(r["weight"] for r in results)
    print("-" * 80)
    print(f"{'TOTAL':<25} {'':<10} {'':<10} {'':<8} {'':<10} {total_weight:<8.2f}")

    print("\nWeights by reporting group:")
    weights_dict = {r["param"]: r["weight"] for r in results}
    for grp, params in REPORT_GROUPS.items():
        total_w = sum(weights_dict.get(p, 0) for p in params)
        details = ", ".join(f"{p}={weights_dict.get(p,0):.1f}" for p in params)
        print(f"  {grp:25s}: {total_w:6.2f}  ({details})")

    # Save to JSON
    output = {
        "schema_version": "4.0-per-parameter-shannon",
        "generated_at": datetime.datetime.utcnow().strftime("%Y-%m-%d"),
        "weight_scheme": "shannon-log-odds-per-parameter",
        "method": {
            "_documentation": [
                "Per-parameter penalty weights using pure Shannon (log-odds) approach.",
                "",
                "For each of the 9 parameters individually:",
                "  - Compute mean Cerny severity in GOOD set vs POOR set",
                "  - weight_i = log(s_poor_i / s_good_i)  [clipped to >= 0]",
                "  - Renormalize to sum = 100",
                "",
                "Cerny severity (Cerny et al. NAR 2026 gkaf1335):",
                "  - severity = 0 if ProSco >= 5 (Preferred)",
                "  - severity = min(1, |Z'|/5) otherwise",
                "",
                "At runtime:",
                "    penalty = Σ penalty_weights[param] × severity[param]",
                "    score   = max(0, min(100, 100 − penalty))",
            ],
            "severity_definition": "Cerny: 0 if ProSco >= 5, else min(1, |Z'|/5)",
            "prosco_preferred_threshold": 5.0,
        },
        "good_reference_set": {
            "n_pairs": n_good,
            "feature_dir": str(args.features_dir),
            "distributions": str(args.distributions),
            "mean_severity": {r["param"]: round(r["mean_good"], 6) for r in results},
        },
        "poor_reference_set": {
            "n_pairs": n_poor,
            "feature_dir": str(args.poor_features_dir),
            "mean_severity": {r["param"]: round(r["mean_poor"], 6) for r in results},
        },
        "raw_logodds_nats": {r["param"]: round(r["raw_logodds"], 4) for r in results},
        "penalty_weights": {r["param"]: round(r["weight"], 4) for r in results},
        "features": FEATURES,
        "report_groups": REPORT_GROUPS,
        "summary_table": [
            {
                "param": r["param"],
                "mean_sev_good": round(r["mean_good"], 4),
                "mean_sev_poor": round(r["mean_poor"], 4),
                "ratio_poor_over_good": round(r["ratio"], 2),
                "raw_logodds": round(r["raw_logodds"], 4),
                "weight": round(r["weight"], 2),
            }
            for r in results
        ],
    }

    with open(args.output, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
