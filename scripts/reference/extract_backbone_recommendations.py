"""Per-residue backbone recommendations over the non-redundant set, with the
suggested correction applied and re-scored.

PARSE flags a low-suiteness residue and names the torsions to correct toward a
target conformer. This asks whether taking that advice works: for every flagged
residue we snap the named torsions to their target and re-run the Richardson
classifier, recording suiteness before and after.

The C++ engine supplies the recommendations and the raw torsions (fast); the
Python classifier does the re-scoring, since only it can classify a suite
vector. The suite is rebuilt here from the per-residue torsions, so the script
also checks its own reconstruction: `suiteness_check` is the classifier's
suiteness for the UNMODIFIED suite, which must match the engine's `suiteness`.
Rows where it does not are dropped and counted.

    python extract_backbone_recommendations.py --workers 8

Writes data/reference/backbone_recommendations.csv (one row per flagged
residue) plus a small JSON summary alongside it.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
ANALYSIS = HERE.parent.parent.parent                 # PARSE-RNA-Analysis/
PROJECTS = ANALYSIS.parent
ENGINE_PY = PROJECTS / "prototyped-pair-finder-main"
PARSE_BIN = PROJECTS / "pair-finder-main-cpp" / "build" / "parse"
CIF_DIR = ENGINE_PY / "data" / "unique_pdbs"
OUT_CSV = ANALYSIS / "data" / "reference" / "backbone_recommendations.csv"
OUT_JSON = ANALYSIS / "data" / "reference" / "backbone_recommendations_summary.json"

sys.path.insert(0, str(ENGINE_PY / "src"))
from pair_finder.scoring.richardson import (  # noqa: E402
    SUITE_ANGLE_NAMES, RichardsonClassifier,
)

# suite = predecessor's (delta, epsilon, zeta) then this residue's
# (alpha, beta, gamma, delta). Mirrors SUITE_ANGLE_NAMES.
assert SUITE_ANGLE_NAMES == ("delta_prev", "epsilon_prev", "zeta_prev",
                             "alpha", "beta", "gamma", "delta")

_CLF: RichardsonClassifier | None = None


def _init():
    global _CLF
    _CLF = RichardsonClassifier()


def _split_res_id(res_id: str):
    """`R-A-20` -> ('R', 20). Insertion codes make the number non-numeric."""
    parts = res_id.split("-")
    if len(parts) < 3:
        return None, None
    try:
        return parts[0], int(parts[-1])
    except ValueError:
        return None, None


def _suite_vector(cur: dict, prev: dict) -> np.ndarray | None:
    vals = [prev.get("delta"), prev.get("epsilon"), prev.get("zeta"),
            cur.get("alpha"), cur.get("beta"), cur.get("gamma"), cur.get("delta")]
    if any(v is None for v in vals):
        return None
    return np.array([float(v) % 360.0 for v in vals])


def _worker(cif: str):
    try:
        raw = subprocess.run([str(PARSE_BIN), cif, "--details", "--no-download"],
                             capture_output=True, text=True, timeout=300)
        if raw.returncode != 0:
            return None, f"parse exit {raw.returncode}"
        doc = json.loads(raw.stdout)
    except Exception as exc:                                    # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}"

    pdb_id = Path(cif).stem.upper()
    recs = doc.get("backbone_residues") or []
    if not recs:
        return [], None

    # index torsions by (chain, seq) so the predecessor is findable
    by_key, by_id = {}, {}
    for t in doc.get("backbone_torsions") or []:
        ch, num = _split_res_id(t["res_id"])
        if ch is not None:
            by_key[(ch, num)] = t
        by_id[t["res_id"]] = t

    rows = []
    for r in recs:
        cur = by_id.get(r["res_id"])
        ch, num = _split_res_id(r["res_id"])
        prev = by_key.get((ch, num - 1)) if ch is not None else None
        if cur is None or prev is None:
            continue
        suite = _suite_vector(cur, prev)
        if suite is None:
            continue

        check = _CLF.classify(suite)
        devs = r.get("deviations") or []
        fixed = suite.copy()
        for d in devs:
            fixed[SUITE_ANGLE_NAMES.index(d["angle"])] = float(d["target"]) % 360.0
        after = _CLF.classify(fixed)

        rows.append({
            "pdb_id": pdb_id,
            "res_id": r["res_id"],
            "base_type": cur.get("base_type"),
            "tier": r.get("tier"),
            "target_conformer": r.get("target_conformer"),
            "suiteness": r.get("suiteness"),
            "suiteness_check": None if check.suiteness is None else round(check.suiteness, 3),
            "conformer_before": check.conformer,
            "suiteness_after": None if after.suiteness is None else round(after.suiteness, 3),
            "conformer_after": after.conformer,
            "n_fired": len(devs),
            "fired_angles": "|".join(d["angle"] for d in devs),
            "worst_prosco": min((d["prosco"] for d in devs if "prosco" in d), default=None),
            "max_abs_gap": max((abs(d["gap"]) for d in devs), default=None),
        })
    return rows, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0, help="first N structures only")
    args = ap.parse_args()

    cifs = sorted(str(p) for p in CIF_DIR.glob("*.cif"))
    if args.limit:
        cifs = cifs[:args.limit]
    print(f"{len(cifs)} structures, {args.workers} workers")

    rows, failed = [], []
    with ProcessPoolExecutor(args.workers, initializer=_init) as ex:
        futs = {ex.submit(_worker, c): c for c in cifs}
        for n, f in enumerate(as_completed(futs), 1):
            got, err = f.result()
            if err:
                failed.append((Path(futs[f]).stem, err))
            else:
                rows.extend(got)
            if n % 250 == 0:
                print(f"  {n}/{len(cifs)}  rows={len(rows)}  failed={len(failed)}")

    df = pd.DataFrame(rows)
    # Reconstruction check: drop rows where our rebuilt suite disagrees with the
    # engine's own suiteness. A mismatch means the predecessor was picked wrong.
    ok = df.suiteness_check.notna() & (
        (df.suiteness - df.suiteness_check).abs() <= 0.01)
    n_bad = int((~ok).sum())
    df = df[ok].drop(columns=["suiteness_check"])

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)

    summary = {
        "n_structures": len(cifs),
        "n_failed": len(failed),
        "n_flagged_residues": int(len(df)),
        "n_dropped_reconstruction_mismatch": n_bad,
        "tier_counts": df.tier.value_counts().to_dict(),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"wrote {OUT_CSV}")
    for pdb, err in failed[:10]:
        print(f"  failed {pdb}: {err}")


if __name__ == "__main__":
    main()
