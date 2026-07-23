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
import os
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
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
DEFAULT_GOOD_SET = REPO_ROOT / "data" / "reference" / "reference_sets" / "reference_set.json"
DEFAULT_POOR_SET = REPO_ROOT / "data" / "reference" / "reference_sets" / "reference_set_poor.json"
# Raw per-structure JSON (only needed for the optional end-to-end Scorer pass) —
# not shipped; set PARSE_JSON_DIR to that directory to use it.
DEFAULT_JSON_DIR = Path(os.environ["PARSE_JSON_DIR"]) if os.environ.get("PARSE_JSON_DIR") else None

# Number of equal-width bins for mutual-information bucket weight derivation.
BUCKET_MI_N_BINS = 20

import _pairfinder_path  # noqa: F401,E402  (adds pair_finder from PARSE_CODE_REPO)
from pair_finder.scoring.canonical import canonicalize_bp, maybe_flip  # noqa: E402
from pair_finder.scoring.hbond_distance_policy import distances_to_score  # noqa: E402
from pair_finder.scoring.issues import REPORT_GROUPS  # noqa: E402

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
# Bucket-level weights (W_PAIRS, W_RESIDUES) via Shannon mutual information
#
# Runs the production Scorer end-to-end on the curated good/poor reference
# PDB sets and computes MI(score_bin; label) for each bucket.  The bucket
# weight is then  w_b = MI_b / Σ MI_*.  Same Shannon-information family as
# the per-parameter penalty weights, generalised from a binary event to a
# continuous score.
# ---------------------------------------------------------------------------

# Module-global for worker processes (set by _init_worker).
_JSON_DIR_FOR_WORKERS: Path | None = None


def _init_worker(json_dir: Path) -> None:
    global _JSON_DIR_FOR_WORKERS
    _JSON_DIR_FOR_WORKERS = json_dir


def _score_one_pdb(args: tuple[str, str]) -> dict | None:
    """Score one PDB end-to-end; return pair_mean + res_mean. Top-level for pickling."""
    pdb_id, bucket = args
    json_dir = _JSON_DIR_FOR_WORKERS
    if json_dir is None or not (json_dir / "pdb_atoms" / f"{pdb_id}.json").exists():
        return None
    try:
        from pair_finder.finder import FinderConfig, PairFinder
        from pair_finder.scoring.scorer import Scorer
        scorer = Scorer.load_default()
        finder = PairFinder(json_dir, FinderConfig(min_score=0.0))
        result = finder.find_pairs(pdb_id, score=False, include_context=True)
        sr = scorer.score_finder_result(result, pdb_id=pdb_id)
    except Exception:
        return None
    return {
        "bucket": bucket,
        "pair_mean": (sum(s.score for s in sr.pair_scores) / len(sr.pair_scores)
                      if sr.pair_scores else None),
        "res_mean": (sum(s.score for s in sr.residue_scores) / len(sr.residue_scores)
                     if sr.residue_scores else None),
    }


def _bin_score(score: float, n_bins: int) -> int:
    b = int(score / (100.0 / n_bins))
    return max(0, min(n_bins - 1, b))


def mutual_information(values_good: list, values_poor: list,
                       n_bins: int = BUCKET_MI_N_BINS) -> tuple:
    """MI(score_bin; label) in bits. Laplace smoothing (eps=0.5/bin)."""
    n_g, n_p = len(values_good), len(values_poor)
    counts_g = [0] * n_bins
    counts_p = [0] * n_bins
    for v in values_good:
        counts_g[_bin_score(v, n_bins)] += 1
    for v in values_poor:
        counts_p[_bin_score(v, n_bins)] += 1
    eps = 0.5
    p_g = [(c + eps) / (n_g + eps * n_bins) for c in counts_g]
    p_p = [(c + eps) / (n_p + eps * n_bins) for c in counts_p]
    n_tot = n_g + n_p
    p_lbl_g, p_lbl_p = n_g / n_tot, n_p / n_tot
    mi = 0.0
    for b in range(n_bins):
        p_bin = p_g[b] * p_lbl_g + p_p[b] * p_lbl_p
        if p_bin <= 0:
            continue
        if p_g[b] * p_lbl_g > 0:
            mi += p_g[b] * p_lbl_g * math.log2(p_g[b] * p_lbl_g / (p_bin * p_lbl_g))
        if p_p[b] * p_lbl_p > 0:
            mi += p_p[b] * p_lbl_p * math.log2(p_p[b] * p_lbl_p / (p_bin * p_lbl_p))
    kl_pg = sum(p_p[b] * math.log2(p_p[b] / p_g[b])
                for b in range(n_bins) if p_p[b] > 0 and p_g[b] > 0)
    kl_gp = sum(p_g[b] * math.log2(p_g[b] / p_p[b])
                for b in range(n_bins) if p_g[b] > 0 and p_p[b] > 0)
    return mi, kl_pg, kl_gp


def compute_bucket_weights(good_set: Path, poor_set: Path, json_dir: Path,
                           workers: int) -> dict:
    good_pdbs = list(json.loads(good_set.read_text())["pdbs"].keys())
    poor_pdbs = list(json.loads(poor_set.read_text())["pdbs"].keys())
    print(f"\n[bucket weights] GOOD={len(good_pdbs)}  POOR={len(poor_pdbs)}  "
          f"workers={workers}  json_dir={json_dir}")

    jobs = [(p, "good") for p in good_pdbs] + [(p, "poor") for p in poor_pdbs]
    t0 = time.time()
    scored: list[dict] = []
    n_err = 0
    with ProcessPoolExecutor(max_workers=workers,
                             initializer=_init_worker,
                             initargs=(json_dir,)) as pool:
        futs = {pool.submit(_score_one_pdb, j): j for j in jobs}
        n_done = 0
        for fut in as_completed(futs):
            n_done += 1
            try:
                r = fut.result(timeout=180)
            except Exception:
                n_err += 1
                continue
            if r is None:
                n_err += 1
                continue
            scored.append(r)
            if n_done % 200 == 0:
                print(f"  [{n_done}/{len(jobs)}] elapsed={time.time()-t0:.0f}s")

    g = [r for r in scored if r["bucket"] == "good"]
    p = [r for r in scored if r["bucket"] == "poor"]
    g_pair = [r["pair_mean"] for r in g if r["pair_mean"] is not None]
    g_res  = [r["res_mean"]  for r in g if r["res_mean"]  is not None]
    p_pair = [r["pair_mean"] for r in p if r["pair_mean"] is not None]
    p_res  = [r["res_mean"]  for r in p if r["res_mean"]  is not None]
    print(f"  Scored {len(scored)} PDBs (errors/missing: {n_err})")
    print(f"  PAIR:    GOOD n={len(g_pair):<5d}  POOR n={len(p_pair)}")
    print(f"  RESIDUE: GOOD n={len(g_res):<5d}  POOR n={len(p_res)}")

    mi_pair, kl_pair_pg, kl_pair_gp = mutual_information(g_pair, p_pair)
    mi_res,  kl_res_pg,  kl_res_gp  = mutual_information(g_res,  p_res)
    total_mi = mi_pair + mi_res
    w_pair = mi_pair / total_mi if total_mi > 0 else 0.5
    w_res  = 1.0 - w_pair

    def _mean(xs):
        return sum(xs) / len(xs) if xs else 0.0

    print(f"\n  MI_pair    = {mi_pair:.4f} bits  → W_PAIRS    = {w_pair:.4f}")
    print(f"  MI_residue = {mi_res:.4f} bits  → W_RESIDUES = {w_res:.4f}")

    return {
        "_documentation": [
            "Bucket weights for overall = W_PAIRS·mean(pair) + W_RESIDUES·mean(residue).",
            "Derived via mutual information MI(score_bin; good/poor label) for each bucket.",
            "  w_b = MI_b / Σ MI_*",
            "Same Shannon-information family as per-parameter penalty weights.",
            "20 equal-width bins on [0,100]; Laplace smoothing eps=0.5/bin.",
        ],
        "W_PAIRS":    round(w_pair, 4),
        "W_RESIDUES": round(w_res,  4),
        "n_bins": BUCKET_MI_N_BINS,
        "mutual_information_bits": {
            "MI_pair":    round(mi_pair, 4),
            "MI_residue": round(mi_res,  4),
        },
        "kl_divergences_bits": {
            "pair":    {"poor||good": round(kl_pair_pg, 4), "good||poor": round(kl_pair_gp, 4)},
            "residue": {"poor||good": round(kl_res_pg,  4), "good||poor": round(kl_res_gp,  4)},
        },
        "supporting_metrics": {
            "pair_bucket":    {"mean_good": round(_mean(g_pair), 2), "mean_poor": round(_mean(p_pair), 2),
                               "n_good": len(g_pair), "n_poor": len(p_pair)},
            "residue_bucket": {"mean_good": round(_mean(g_res),  2), "mean_poor": round(_mean(p_res),  2),
                               "n_good": len(g_res),  "n_poor": len(p_res)},
        },
        "reference_sets": {
            "good": str(good_set), "poor": str(poor_set),
            "json_dir": str(json_dir), "load_errors_or_missing": n_err,
        },
    }


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
    parser.add_argument(
        "--good-set",
        type=Path,
        default=DEFAULT_GOOD_SET,
        help="reference_set.json (good PDBs) for bucket weight derivation",
    )
    parser.add_argument(
        "--poor-set",
        type=Path,
        default=DEFAULT_POOR_SET,
        help="reference_set_poor.json for bucket weight derivation",
    )
    parser.add_argument(
        "--json-dir",
        type=Path,
        default=DEFAULT_JSON_DIR,
        help="raw per-structure JSON dir (PARSE_JSON_DIR) for the optional end-to-end Scorer pass",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=8,
        help="Parallel workers for the end-to-end Scorer pass (default: 8)",
    )
    parser.add_argument(
        "--skip-bucket-weights",
        action="store_true",
        help="Skip the bucket weight pass (per-parameter weights only, ~10x faster)",
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

    # ------------------------------------------------------------------
    # Bucket-level weights (W_PAIRS, W_RESIDUES)
    # Runs the production Scorer end-to-end on good/poor PDB sets.
    # The per-parameter weights written above are loaded by each worker
    # via Scorer.load_default(), so both derivations are consistent.
    # ------------------------------------------------------------------
    if args.skip_bucket_weights:
        print("\n[bucket weights] Skipped (--skip-bucket-weights).")
        return

    bucket_report = compute_bucket_weights(
        args.good_set, args.poor_set, args.json_dir, args.workers,
    )
    output["bucket_weights"] = bucket_report
    with open(args.output, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nRewrote {args.output} with bucket_weights section.")
    print(f"  W_PAIRS    = {bucket_report['W_PAIRS']}")
    print(f"  W_RESIDUES = {bucket_report['W_RESIDUES']}")


if __name__ == "__main__":
    main()
