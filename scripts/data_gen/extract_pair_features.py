"""Extract per-pair, per-hbond, per-residue features from the reference set.

Loops over every PDB in data/reference/reference_sets/reference_set.json, runs
the PARSE C++ engine (`parse --details`), and writes three CSV tables under
data/reference/high_quality_features/:

    pairs.csv     — one row per base-base pair (Tsukuba params, plane angle,
                        template RMSD, lw_class, etc.)
    hbonds.csv    — one row per H-bond inside a base-base pair (distance,
                        donor / acceptor angles, slot alignment)
    torsions.csv  — one row per residue-in-a-pair (alpha-zeta + Richardson
                        suite), tagged with the pair's bp_type / lw_class

Non-base-base or is_ambiguous pairs are excluded (they would muddy the
per-(bp_type, edge_type) reference distributions). The engine downloads each
mmCIF itself, so no local structure files are needed — see _engine.py.
"""

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

import _engine
from _engine import is_scorable, res_name, run_parse

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_REF = REPO_ROOT / "data" / "reference" / "reference_sets" / "reference_set.json"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "data" / "reference" / "high_quality_features"
PARAMS = ["shear", "stretch", "stagger", "buckle", "propeller", "opening"]


def extract_one(pdb_id: str, d: dict) -> tuple[list, list, list]:
    pairs, hbonds, torsions = [], [], []
    bt_by_res = {r["res_id"]: r for r in d.get("backbone_torsions", [])}
    idx = 0
    for p in d["pairs"]:
        if not is_scorable(p):
            continue
        bp_type = p["bp_type"]
        lw = p["classification"]["lw_class"]
        rb = p["rigid_body"]
        n_base = sum(1 for c in p["hbond_categories"] if c == "base-base")
        pairs.append({
            "pdb_id": pdb_id, "pair_idx": idx,
            "res_id1": p["res_id1"], "res_id2": p["res_id2"],
            "res_name1": res_name(p["res_id1"]), "res_name2": res_name(p["res_id2"]),
            "bp_type": bp_type, "lw_class": lw,
            **{k: rb[k] for k in PARAMS},
            "plane_angle": p["geometry"]["plane_angle"],
            "template_rmsd": p.get("template_rmsd"),
            "num_hbonds": p["num_hbonds"], "num_base_hbonds": n_base,
            "quality_score": p["quality_score"],
        })
        for hb in p.get("hbonds", []):
            hbonds.append({
                "pdb_id": pdb_id, "pair_idx": idx, "bp_type": bp_type, "lw_class": lw,
                "donor_atom": hb["donor_atom"], "acceptor_atom": hb["acceptor_atom"],
                "distance": hb["distance"], "donor_angle": hb["donor_angle"],
                "acceptor_angle": hb["acceptor_angle"], "slot_alignment": hb["slot_alignment"],
            })
        for pos, rid in enumerate([p["res_id1"], p["res_id2"]], 1):
            bt = bt_by_res.get(rid)
            if bt is None:
                continue
            torsions.append({
                "pdb_id": pdb_id, "pair_idx": idx, "res_id": rid,
                "base_type": bt["base_type"], "bp_type": bp_type, "lw_class": lw,
                "residue_position": pos,
                "alpha": bt["alpha"], "beta": bt["beta"], "gamma": bt["gamma"],
                "delta": bt["delta"], "epsilon": bt["epsilon"], "zeta": bt["zeta"],
                "conformer": bt["conformer"], "suiteness": bt["suiteness"],
                "suite_bin": bt["suite_bin"], "is_outlier": bt["is_outlier"],
            })
        idx += 1
    return pairs, hbonds, torsions


_ENGINE: str | None = None


def _init_worker(engine: str) -> None:
    global _ENGINE
    _ENGINE = engine


def _worker(pdb_id: str):
    try:
        return pdb_id, *extract_one(pdb_id, run_parse(pdb_id, _ENGINE)), None
    except Exception as e:  # noqa: BLE001
        return pdb_id, [], [], [], f"{type(e).__name__}: {e}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference-set", type=Path, default=DEFAULT_REF)
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    ap.add_argument("--engine", default=None, help="path to `parse` (else PARSE_ENGINE / PATH)")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    engine = _engine.find_engine(args.engine)
    with open(args.reference_set) as f:
        pdbs = list(json.load(f)["pdbs"])
    if args.limit:
        pdbs = pdbs[:args.limit]

    print(f"Engine:          {engine}")
    print(f"Reference set:   {args.reference_set}")
    print(f"PDBs to extract: {len(pdbs)}")
    print(f"Output:          {args.output_dir}/")

    all_pairs, all_hbonds, all_torsions, failed = [], [], [], []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers,
                             initializer=_init_worker, initargs=(engine,)) as ex:
        for i, (pdb_id, pr, hb, to, err) in enumerate(ex.map(_worker, pdbs), 1):
            if err is not None:
                failed.append((pdb_id, err))
            else:
                all_pairs += pr; all_hbonds += hb; all_torsions += to
            if i % 50 == 0 or i == len(pdbs):
                rate = i / (time.time() - t0)
                print(f"  [{i}/{len(pdbs)}] pairs={len(all_pairs):,} hbonds={len(all_hbonds):,} "
                      f"torsions={len(all_torsions):,} failed={len(failed)} | {rate:.1f} pdb/s", flush=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(all_pairs).to_csv(args.output_dir / "pairs.csv", index=False)
    pd.DataFrame(all_hbonds).to_csv(args.output_dir / "hbonds.csv", index=False)
    pd.DataFrame(all_torsions).to_csv(args.output_dir / "torsions.csv", index=False)
    summary = {
        "n_pdbs_attempted": len(pdbs), "n_pdbs_failed": len(failed),
        "n_pairs": len(all_pairs), "n_hbonds": len(all_hbonds),
        "n_torsion_observations": len(all_torsions),
        "elapsed_seconds": round(time.time() - t0, 1),
        "failed_pdbs": [{"pdb_id": p, "error": e} for p, e in failed],
    }
    with open(args.output_dir / "extraction_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nWrote feature tables to {args.output_dir}/")
    print(f"  pairs.csv    : {len(all_pairs):,} rows")
    print(f"  hbonds.csv   : {len(all_hbonds):,} rows")
    print(f"  torsions.csv : {len(all_torsions):,} rows")
    if failed:
        print(f"  failed PDBs  : {len(failed)} (see extraction_summary.json)")


if __name__ == "__main__":
    main()
