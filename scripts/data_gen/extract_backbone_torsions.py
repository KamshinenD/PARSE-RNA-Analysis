"""Extract per-residue backbone suites for EVERY classifiable residue (paired +
unpaired) in the high-quality reference set.

Unlike extract_pair_features.py (one torsion row per residue-in-a-pair), this
walks every residue with a valid 7D suite in each structure, so the
per-conformer torsion distributions (build_backbone_distributions.py) represent
loops/junctions — where backbone errors concentrate — not just helical context.

Runs the PARSE C++ engine (`parse --details`) and reads its top-level
`backbone_torsions` array. The engine downloads each mmCIF itself (see
_engine.py).

Output: data/reference/high_quality_features/torsions_all.csv
    pdb_id, res_id, base_type, alpha..zeta, conformer, suiteness, suite_bin, is_outlier
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from pathlib import Path

import pandas as pd

import _engine
from _engine import run_parse

REPO = Path(__file__).resolve().parent.parent.parent
DEFAULT_REF = REPO / "data" / "reference" / "reference_sets" / "reference_set.json"
DEFAULT_OUT = REPO / "data" / "reference" / "high_quality_features" / "torsions_all.csv"
COLS = ["res_id", "base_type", "alpha", "beta", "gamma", "delta", "epsilon", "zeta",
        "conformer", "suiteness", "suite_bin", "is_outlier"]

_ENGINE: str | None = None


def _init_worker(engine: str) -> None:
    global _ENGINE
    _ENGINE = engine


def _worker(pdb_id: str):
    try:
        d = run_parse(pdb_id, _ENGINE)
        rows = [{"pdb_id": pdb_id, **{c: r[c] for c in COLS}} for r in d.get("backbone_torsions", [])]
        return pdb_id, rows, None
    except Exception as e:  # noqa: BLE001
        return pdb_id, [], f"{type(e).__name__}: {e}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference-set", type=Path, default=DEFAULT_REF)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--engine", default=None, help="path to `parse` (else PARSE_ENGINE / PATH)")
    ap.add_argument("--workers", type=int, default=0, help="0 = cpu_count")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    engine = _engine.find_engine(args.engine)
    with open(args.reference_set) as f:
        pdbs = list(json.load(f)["pdbs"])
    if args.limit:
        pdbs = pdbs[:args.limit]
    print(f"Engine: {engine}\nReference set: {len(pdbs)} structures")

    all_rows: list[dict] = []
    n_ok = n_fail = 0
    with mp.Pool(args.workers or mp.cpu_count(),
                 initializer=_init_worker, initargs=(engine,)) as pool:
        for i, (pdb_id, rows, err) in enumerate(pool.imap_unordered(_worker, pdbs), 1):
            if err:
                n_fail += 1
            else:
                all_rows.extend(rows); n_ok += 1
            if i % 50 == 0 or i == len(pdbs):
                print(f"  {i}/{len(pdbs)}  ok={n_ok} fail={n_fail}  rows={len(all_rows):,}", flush=True)

    df = pd.DataFrame(all_rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"\nWrote {args.output}")
    print(f"  {len(df):,} residues  |  {df.pdb_id.nunique()} structures  |  "
          f"{df.conformer.nunique()} conformers  |  in-cluster {int((~df['is_outlier']).sum()):,}")


if __name__ == "__main__":
    main()
