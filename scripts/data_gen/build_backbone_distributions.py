"""Build per-(conformer, suite-torsion) ProSco + Z' tables from the GOOD reference.

Backbone analogue of build_prosco_distributions.py / build_z_tables.py. The
scoring *cell* is the Richardson conformer (the backbone twin of bp_type×lw_class):
a torsion is only unimodal — and therefore ProSco/Z'-scorable — *within* a
conformer. Pooling conformers makes each torsion multimodal, so there is NO
_ANY fallback here; conformers with too few observations are simply omitted and
the runtime falls back to the Richardson fixed width.

Two backbone-specific details vs the pair pipeline:

  * Suite reconstruction. torsions.csv stores per-residue alpha..zeta; the
    7D suite needs the predecessor's delta/epsilon/zeta. We dedup to one row per
    residue, join each residue's backbone predecessor, and assemble
    [delta_prev, epsilon_prev, zeta_prev, alpha, beta, gamma, delta] in 0-360.

  * Circular unwrap. Torsions wrap at 0/360, so before any linear KDE/quantile
    we shift each value onto the branch nearest its conformer center:
        unwrapped = center + signed_gap(value, center)
    Values then sit contiguously in (center-180, center+180]. The center is
    stored per cell as ``unwrap_center`` so the runtime reproduces it exactly.

Outputs:
    data/reference/scoring_tables/backbone_prosco_distributions.json   {"backbone": {angle: {conf: cell}}}
    data/reference/scoring_tables/backbone_z_tables.json               {"backbone": {angle: {conf: cell}}}
"""

from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "scripts" / "data_gen"))
from _vendor.helpers import parse_res_seq  # noqa: E402
from _vendor.richardson import (  # noqa: E402
    OUTLIER_CONFORMER, SUITE_ANGLE_NAMES, RichardsonClassifier, signed_angle_gap,
)
from build_prosco_distributions import build_prosco_cell  # noqa: E402

REF = REPO / "data" / "reference"
PHI_INV_09 = 1.2815515655446004
MIN_OBS = 10

# The three suite positions taken from the residue itself vs its predecessor.
# suite = [delta_prev, epsilon_prev, zeta_prev, alpha, beta, gamma, delta]
_CUR_COLS = {"alpha": 3, "beta": 4, "gamma": 5, "delta": 6}
_PREV_COLS = {"delta": 0, "epsilon": 1, "zeta": 2}


def reconstruct_suites(df: pd.DataFrame) -> pd.DataFrame:
    """One row per residue with its 7D suite (0-360), conformer, suite_bin, pdb.

    Keeps only non-outlier residues with a real conformer — the reference
    distribution of where each torsion sits *when the residue is in* that cluster.
    """
    res = df.drop_duplicates(subset=["pdb_id", "res_id"]).copy()

    rows = []
    for pdb_id, g in res.groupby("pdb_id"):
        by_id = {r.res_id: r for r in g.itertuples(index=False)}
        # order residues within each chain by (seq_num, icode)
        ordered: dict[str, list] = {}
        for r in g.itertuples(index=False):
            p = parse_res_seq(r.res_id)
            if p is None:
                continue
            ordered.setdefault(p[0], []).append((p[1], p[2], r.res_id))
        pred: dict[str, str] = {}
        for chain_res in ordered.values():
            chain_res.sort(key=lambda t: (t[0], t[1]))
            for i in range(1, len(chain_res)):
                pn, _, pid = chain_res[i - 1]
                cn, _, cid = chain_res[i]
                if pn == cn or pn == cn - 1:
                    pred[cid] = pid

        for r in g.itertuples(index=False):
            if r.is_outlier or r.conformer in (None, OUTLIER_CONFORMER) or pd.isna(r.conformer):
                continue
            pid = pred.get(r.res_id)
            if pid is None:
                continue
            pr = by_id[pid]
            vals = [pr.delta, pr.epsilon, pr.zeta, r.alpha, r.beta, r.gamma, r.delta]
            if any(v is None or (isinstance(v, float) and np.isnan(v)) for v in vals):
                continue
            suite = [float(v) % 360.0 for v in vals]
            rows.append((pdb_id, r.conformer, r.suite_bin, *suite))

    cols = ["pdb_id", "conformer", "suite_bin", *SUITE_ANGLE_NAMES]
    return pd.DataFrame(rows, columns=cols)


def per_pdb_weights(pdb_ids: pd.Series) -> np.ndarray:
    counts = pdb_ids.value_counts()
    return np.array([1.0 / counts[p] for p in pdb_ids])


def weighted_quantile(values, weights, q) -> float:
    order = np.argsort(values)
    v, w = values[order], weights[order]
    cum = np.cumsum(w)
    return float(np.interp(q * cum[-1], cum - w / 2.0, v))


def build_z_cell(values: np.ndarray, weights: np.ndarray, center: float) -> dict | None:
    if len(values) < MIN_OBS or np.ptp(values) < 1e-9:
        return None
    M = weighted_quantile(values, weights, 0.5)
    D1 = weighted_quantile(values, weights, 0.1)
    D9 = weighted_quantile(values, weights, 0.9)
    return {
        "M": round(M, 4),
        "sl": round(max((M - D1) / PHI_INV_09, 1e-9), 4),
        "su": round(max((D9 - M) / PHI_INV_09, 1e-9), 4),
        "n": int(len(values)),
        "unwrap_center": round(center, 4),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--torsions", type=Path,
                    default=REF / "high_quality_features" / "torsions_all.csv",
                    help="per-residue torsions CSV (all residues, paired+unpaired)")
    ap.add_argument("--prosco-out", type=Path, default=REF / "scoring_tables" / "backbone_prosco_distributions.json")
    ap.add_argument("--z-out", type=Path, default=REF / "scoring_tables" / "backbone_z_tables.json")
    args = ap.parse_args()

    print(f"Reading {args.torsions}")
    df = pd.read_csv(args.torsions)
    print(f"  {len(df):,} residue rows / {df.pdb_id.nunique()} structures")

    print("Reconstructing 7D suites (dedup + predecessor join) ...")
    suites = reconstruct_suites(df)
    print(f"  {len(suites):,} classified suites across {suites.conformer.nunique()} conformers")

    clf = RichardsonClassifier()
    prosco: dict[str, dict[str, dict]] = {a: {} for a in SUITE_ANGLE_NAMES}
    ztab: dict[str, dict[str, dict]] = {a: {} for a in SUITE_ANGLE_NAMES}

    n_cells = 0
    for conformer, g in suites.groupby("conformer"):
        center = clf.center_for(conformer)
        if center is None:
            continue
        w = per_pdb_weights(g["pdb_id"])
        pdb_series = g["pdb_id"].reset_index(drop=True)
        for k, angle in enumerate(SUITE_ANGLE_NAMES):
            c = float(center[k])
            # circular unwrap onto the branch nearest the conformer center
            unwrapped = np.array([c + signed_angle_gap(v, c) for v in g[angle].to_numpy(float)])
            pcell = build_prosco_cell(unwrapped, w, pdb_series)
            if pcell["prosco_per_bin"] is not None:
                pcell["unwrap_center"] = round(c, 4)
                prosco[angle][conformer] = pcell
            zcell = build_z_cell(unwrapped, w, c)
            if zcell is not None:
                ztab[angle][conformer] = zcell
        n_cells += 1

    prosco_out = {
        "schema_version": "1.0",
        "generated_at": datetime.date.today().isoformat(),
        "method": {
            "cell": "Richardson conformer (no _ANY fallback: pooling conformers is multimodal)",
            "circular": "values unwrapped as center + signed_gap(value, center) before KDE; "
                        "runtime must unwrap the query the same way using unwrap_center",
            "prosco": "same KDE/Nyquist ProSco as build_prosco_distributions.py",
        },
        "prosco_distributions": {"backbone": prosco},
    }
    z_out = {
        "schema_version": "1.0",
        "generated_at": datetime.date.today().isoformat(),
        "method": {"M": "weighted median (unwrapped)", "sl": "(M-D1)/1.2816",
                   "su": "(D9-M)/1.2816", "circular": "unwrap via unwrap_center"},
        "z_tables": {"backbone": ztab},
    }

    args.prosco_out.write_text(json.dumps(prosco_out, ensure_ascii=False))
    args.z_out.write_text(json.dumps(z_out, indent=2))
    print(f"\nWrote {args.prosco_out.name} and {args.z_out.name}")
    print(f"  conformer cells: {n_cells}")
    for angle in SUITE_ANGLE_NAMES:
        print(f"  {angle:14s}: prosco {len(prosco[angle]):2d} conf  |  z' {len(ztab[angle]):2d} conf")


if __name__ == "__main__":
    main()
