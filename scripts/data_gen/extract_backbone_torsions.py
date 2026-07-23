"""Extract per-residue backbone suites for EVERY residue (paired + unpaired) in
the high-quality reference set.

Unlike extract_pair_features.py (which emits one torsion row per residue-in-a-
pair, so loops are missing), this walks every classifiable residue in each
structure. That makes the per-conformer torsion distributions
(build_backbone_distributions.py) represent loops/junctions — where backbone
errors actually concentrate — not just helical, paired context.

Output: data/reference/high_quality_features/torsions_all.csv
    pdb_id, res_id, base_type, alpha..zeta, conformer, suiteness, suite_bin, is_outlier
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
import _pairfinder_path  # noqa: F401,E402  (adds pair_finder from PARSE_CODE_REPO)
from pair_finder.core.residue import normalize_base_type  # noqa: E402
from pair_finder.finder.cache import PairCache  # noqa: E402
from pair_finder.scoring.features import (  # noqa: E402
    build_predecessor_index, build_suite_for,
)
from pair_finder.scoring.richardson import RichardsonClassifier  # noqa: E402
from scorer_parameters import compute_all_torsions  # noqa: E402

DEFAULT_REF = REPO / "data" / "reference" / "reference_sets" / "reference_set.json"
# Raw per-structure JSON (pair-finder input) — not shipped in this package; set
# PARSE_JSON_DIR to that directory (e.g. on the external drive) to run extraction.
DEFAULT_JSON_DIR = Path(os.environ["PARSE_JSON_DIR"]) if os.environ.get("PARSE_JSON_DIR") else None
DEFAULT_OUT = REPO / "data" / "reference" / "high_quality_features" / "torsions_all.csv"

_W_JSON_DIR: Path | None = None
_W_CLS: RichardsonClassifier | None = None


def _init_worker(json_dir: str) -> None:
    global _W_JSON_DIR, _W_CLS
    _W_JSON_DIR = Path(json_dir)
    _W_CLS = RichardsonClassifier()


def _extract_worker(pdb_id: str):
    """Return (pdb_id, rows, err) — one row per classifiable residue.

    Loads ONLY the atoms (PairCache._load_residues) — no pair-finding, which is
    the expensive part on large structures. Torsions/suiteness are identical
    to the finder path (differential-tested), just far faster on ribosomes.
    """
    rows: list[dict] = []
    try:
        cache = PairCache(pdb_id, _W_JSON_DIR)
        cache._load_residues()
        residues = cache.residues
        all_torsions = compute_all_torsions(residues)
        predecessor = build_predecessor_index(all_torsions)
        for res_id, t in all_torsions.items():
            suite = build_suite_for(res_id, all_torsions, predecessor)
            if suite is None:
                continue
            sr = _W_CLS.classify(suite)
            residue = residues.get(res_id)
            base = normalize_base_type(residue.base_type) if residue else None
            rows.append({
                "pdb_id": pdb_id,
                "res_id": res_id,
                "base_type": base,
                "alpha": t.alpha, "beta": t.beta, "gamma": t.gamma,
                "delta": t.delta, "epsilon": t.epsilon, "zeta": t.zeta,
                "conformer": sr.conformer,
                "suiteness": sr.suiteness,
                "suite_bin": sr.bin,
                "is_outlier": sr.is_outlier,
            })
    except Exception as e:  # noqa: BLE001
        return pdb_id, rows, f"{type(e).__name__}: {e}"
    return pdb_id, rows, None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference-set", type=Path, default=DEFAULT_REF)
    ap.add_argument("--json-dir", type=Path, default=DEFAULT_JSON_DIR)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--workers", type=int, default=0, help="0 = cpu_count")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    with open(args.reference_set) as f:
        ref = json.load(f)
    pdbs = [p for p, e in ref["pdbs"].items() if e.get("json_ready")]
    if args.limit:
        pdbs = pdbs[:args.limit]
    print(f"Reference set: {len(pdbs)} json-ready structures")

    import multiprocessing as mp
    n_workers = args.workers or mp.cpu_count()
    all_rows: list[dict] = []
    n_ok = n_fail = 0
    with mp.Pool(n_workers, initializer=_init_worker,
                 initargs=(str(args.json_dir),)) as pool:
        for i, (pdb_id, rows, err) in enumerate(
                pool.imap_unordered(_extract_worker, pdbs), 1):
            if err:
                n_fail += 1
            else:
                all_rows.extend(rows)
                n_ok += 1
            if i % 50 == 0 or i == len(pdbs):
                print(f"  {i}/{len(pdbs)}  ok={n_ok} fail={n_fail} "
                      f" rows={len(all_rows):,}", flush=True)

    df = pd.DataFrame(all_rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    n_out = int((~df["is_outlier"]).sum())
    print(f"\nWrote {args.output}")
    print(f"  {len(df):,} residues  |  {df.pdb_id.nunique()} structures  |  "
          f"{df.conformer.nunique()} conformers  |  in-cluster {n_out:,}")


if __name__ == "__main__":
    main()
