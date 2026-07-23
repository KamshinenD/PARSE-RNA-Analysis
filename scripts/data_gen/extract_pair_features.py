"""Extract per-pair, per-hbond, per-residue features from the reference set.

Loops over every JSON-ready PDB in data/reference/reference_sets/reference_set.json,
runs pair-finder, computes the geometric / H-bond / torsion features, and
writes them as three CSV tables under data/reference/high_quality_features/:

    pairs.csv     — one row per base-base pair (Tsukuba params, plane angle,
                        template RMSD, lw_class, etc.)
    hbonds.csv    — one row per H-bond inside a base-base pair (distance,
                        donor / acceptor angles, slot alignment)
    torsions.csv  — one row per residue-in-a-pair (alpha-zeta),
                        tagged with the pair's bp_type / lw_class for context

Pairs that are non-base-base (sugar-edge fallbacks, etc.) or flagged as
``is_ambiguous`` are excluded — they would muddy the per-(bp_type,
edge_type) reference distributions.

Modified bases are normalized to their parent base type (PSU → U, 1MA → A,
OMG → G, etc.) for the bp_type field; raw names are preserved alongside.
"""

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

# Package-relative paths (script is at scripts/data_gen/, package root is two up)
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_REF = REPO_ROOT / "data" / "reference" / "reference_sets" / "reference_set.json"
# Raw per-structure JSON (pair-finder input) — not shipped in this package; set
# PARSE_JSON_DIR to that directory (e.g. on the external drive) to run extraction.
DEFAULT_JSON_DIR = Path(os.environ["PARSE_JSON_DIR"]) if os.environ.get("PARSE_JSON_DIR") else None
DEFAULT_OUTPUT_DIR = REPO_ROOT / "data" / "reference" / "high_quality_features"

# Bring src/ onto sys.path so we can import pair_finder + scorer_parameters
# even when this script is launched outside an installed venv.
import _pairfinder_path  # noqa: F401,E402  (adds pair_finder from PARSE_CODE_REPO)

from pair_finder.finder import PairFinder  # noqa: E402
from pair_finder.hbond.finder import HBondFinder  # noqa: E402
from pair_finder.core.residue import normalize_base_type  # noqa: E402
from pair_finder.scoring.richardson import RichardsonClassifier  # noqa: E402
from scorer_parameters import compute_pair_parameters, compute_all_torsions  # noqa: E402
from scorer_parameters.hbond_angles import compute_hbond_angles  # noqa: E402


def _build_predecessor_index(all_torsions: dict) -> dict:
    """Return a (chain, seq_num) → res_id lookup over residues with torsions.

    Used to find each residue's predecessor in the chain so we can build
    the 7D Richardson suite (which needs delta_prev, epsilon_prev, zeta_prev).
    """
    lookup = {}
    for res_id in all_torsions:
        parts = res_id.split("-")
        if len(parts) < 3:
            continue
        try:
            chain = parts[0]
            num = int(parts[2])
        except (ValueError, IndexError):
            continue
        lookup[(chain, num)] = res_id
    return lookup


def _build_suite_for(res_id: str, all_torsions: dict, predecessor_index: dict):
    """Return the 7D suite array for ``res_id`` or None if unavailable.

    Suite ordering (Richardson 2008):
        [delta_prev, epsilon_prev, zeta_prev, alpha, beta, gamma, delta]
    """
    parts = res_id.split("-")
    if len(parts) < 3:
        return None
    try:
        chain, num = parts[0], int(parts[2])
    except ValueError:
        return None

    pred_id = predecessor_index.get((chain, num - 1))
    if pred_id is None:
        return None

    cur = all_torsions.get(res_id)
    pred = all_torsions.get(pred_id)
    if cur is None or pred is None:
        return None

    needed = (pred.delta, pred.epsilon, pred.zeta,
              cur.alpha, cur.beta, cur.gamma, cur.delta)
    if any(v is None for v in needed):
        return None

    import numpy as np  # local to keep import light when unused
    return np.array([v % 360 for v in needed], dtype=float)


def extract_one_pdb(pdb_id: str, finder: PairFinder, hbond_finder: HBondFinder,
                    classifier: RichardsonClassifier):
    """Return three lists of dict records for one PDB:
        pair_records, hbond_records, torsion_records.

    Skips non-base-base pairs and ambiguous LW classifications.
    Each torsion record includes a Richardson suite classification
    (conformer name + suiteness) and a chi class (anti/syn).
    """
    pair_records: list[dict] = []
    hbond_records: list[dict] = []
    torsion_records: list[dict] = []

    result = finder.find_pairs(pdb_id, include_context=True)
    if not result.pairs:
        return pair_records, hbond_records, torsion_records

    # Compute torsions once for all residues in the structure (also used as
    # source for predecessor delta/epsilon/zeta when building the 7D suite).
    all_torsions = compute_all_torsions(result.atoms)
    predecessor_index = _build_predecessor_index(all_torsions)

    for pair_idx, pair in enumerate(result.pairs):
        # Reference set should reflect canonical RNA base pairs only.
        if pair.pair_category != "base-base":
            continue
        if pair.is_ambiguous:
            continue
        if pair.lw_class is None or "|" in pair.lw_class:
            continue

        # Tsukuba parameters from the pair's reference frames.
        params = compute_pair_parameters(pair.frame1, pair.frame2)

        # plane_angle kept for informational output even though it's not
        # scored (its check was dropped — buckle + stagger cover the same
        # signal). d_v dropped entirely since nothing references it.
        v = pair.validation
        plane_angle = float(v.plane_angle) if v is not None else None

        # Normalize modified bases (PSU → U, 1MA → A, ...) for grouping.
        nb1 = normalize_base_type(pair.res_name1)
        nb2 = normalize_base_type(pair.res_name2)
        bp_type = f"{nb1}-{nb2}"

        pair_records.append({
            "pdb_id": pdb_id,
            "pair_idx": pair_idx,
            "res_id1": pair.res_id1,
            "res_id2": pair.res_id2,
            "res_name1": pair.res_name1,
            "res_name2": pair.res_name2,
            "bp_type": bp_type,
            "lw_class": pair.lw_class,
            "shear": params.shear,
            "stretch": params.stretch,
            "stagger": params.stagger,
            "buckle": params.buckle,
            "propeller": params.propeller,
            "opening": params.opening,
            "plane_angle": plane_angle,
            "template_rmsd": pair.template_rmsd,
            "num_hbonds": pair.num_hbonds,
            "num_base_hbonds": pair.num_base_hbonds,
            "quality_score": pair.quality_score,
        })

        # Per-H-bond features (re-run HBondFinder; pair-finder doesn't preserve
        # the per-bond list on the candidate, only counts).
        res1 = result.atoms.get(pair.res_id1)
        res2 = result.atoms.get(pair.res_id2)
        if res1 is not None and res2 is not None:
            hbonds = hbond_finder.find_between(res1, res2)
            for hb in hbonds:
                # Locate which residue carries the donor (for the angle helper).
                if hb.donor_res_id == pair.res_id1:
                    donor_res, acceptor_res = res1, res2
                else:
                    donor_res, acceptor_res = res2, res1
                donor_angle, acceptor_angle = compute_hbond_angles(
                    hb, donor_res, acceptor_res
                )
                hbond_records.append({
                    "pdb_id": pdb_id,
                    "pair_idx": pair_idx,
                    "bp_type": bp_type,
                    "lw_class": pair.lw_class,
                    "donor_atom": hb.donor_atom,
                    "acceptor_atom": hb.acceptor_atom,
                    "distance": hb.distance,
                    "donor_angle": donor_angle,
                    "acceptor_angle": acceptor_angle,
                    "slot_alignment": hb.alignment_score,
                })

        # Per-residue torsions, tagged with this pair's context.
        for pos, res_id in enumerate([pair.res_id1, pair.res_id2], 1):
            t = all_torsions.get(res_id)
            if t is None:
                continue
            base_type = nb1 if pos == 1 else nb2

            # Richardson suite classification (needs predecessor torsions).
            suite = _build_suite_for(res_id, all_torsions, predecessor_index)
            if suite is not None:
                sr = classifier.classify(suite)
                conformer = sr.conformer
                suiteness = sr.suiteness
                suite_bin = sr.bin
                is_outlier = sr.is_outlier
            else:
                # Predecessor missing (chain start or unparseable res_id).
                conformer = None
                suiteness = None
                suite_bin = None
                is_outlier = None

            torsion_records.append({
                "pdb_id": pdb_id,
                "pair_idx": pair_idx,
                "res_id": res_id,
                "base_type": base_type,
                "bp_type": bp_type,
                "lw_class": pair.lw_class,
                "residue_position": pos,
                "alpha": t.alpha,
                "beta": t.beta,
                "gamma": t.gamma,
                "delta": t.delta,
                "epsilon": t.epsilon,
                "zeta": t.zeta,
                "conformer": conformer,
                "suiteness": suiteness,
                "suite_bin": suite_bin,
                "is_outlier": is_outlier,
            })

    return pair_records, hbond_records, torsion_records


# Per-worker finder objects (one set per process; not picklable, so built in
# an initializer rather than passed in). 'spawn' re-imports this module in each
# child, so the module-level imports above are available.
_W_FINDER = _W_HB = _W_CLS = None


def _init_worker(json_dir: str) -> None:
    global _W_FINDER, _W_HB, _W_CLS
    _W_FINDER = PairFinder(json_dir=Path(json_dir))
    _W_HB = HBondFinder()
    _W_CLS = RichardsonClassifier()


def _extract_worker(pdb_id: str):
    """Extract one PDB in a worker. Returns (pdb_id, pairs, hbonds, torsions, err)."""
    try:
        p, h, t = extract_one_pdb(pdb_id, _W_FINDER, _W_HB, _W_CLS)
        return pdb_id, p, h, t, None
    except Exception as e:  # noqa: BLE001
        return pdb_id, [], [], [], f"{type(e).__name__}: {e}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference-set", type=Path, default=DEFAULT_REF)
    ap.add_argument("--json-dir", type=Path, default=DEFAULT_JSON_DIR)
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    ap.add_argument("--workers", type=int, default=6,
                    help="Parallel worker processes (default 6; 1 = sequential)")
    ap.add_argument("--limit", type=int, default=None,
                    help="Process only the first N JSON-ready PDBs (smoke testing)")
    args = ap.parse_args()

    with open(args.reference_set) as f:
        ref = json.load(f)

    pdbs = [p for p, e in ref["pdbs"].items() if e["json_ready"]]
    if args.limit:
        pdbs = pdbs[:args.limit]

    print(f"Reference set:      {args.reference_set}")
    print(f"JSON dir:           {args.json_dir}")
    print(f"PDBs to extract:    {len(pdbs)}")
    print(f"Output:             {args.output_dir}/")
    print(f"Workers:            {args.workers}")

    all_pairs: list[dict] = []
    all_hbonds: list[dict] = []
    all_torsions: list[dict] = []
    failed: list[tuple[str, str]] = []

    # Workers return their records; the parent concatenates and writes the
    # CSV files once at the end, so there is no shared-state race.
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(max_workers=args.workers,
                             initializer=_init_worker,
                             initargs=(str(args.json_dir),)) as ex:
        for pdb_id, p, h, t, err in ex.map(_extract_worker, pdbs):
            done += 1
            if err is not None:
                failed.append((pdb_id, err))
            else:
                all_pairs.extend(p)
                all_hbonds.extend(h)
                all_torsions.extend(t)
            if done % 50 == 0 or done == len(pdbs):
                elapsed = time.time() - t0
                rate = done / elapsed if elapsed > 0 else 0
                eta = (len(pdbs) - done) / rate if rate > 0 else 0
                print(f"  [{done}/{len(pdbs)}] pairs={len(all_pairs):,} "
                      f"hbonds={len(all_hbonds):,} torsions={len(all_torsions):,} "
                      f"failed={len(failed)} | {rate:.1f} pdb/s, ETA {eta:.0f}s",
                      flush=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pairs_df = pd.DataFrame(all_pairs)
    hbonds_df = pd.DataFrame(all_hbonds)
    torsions_df = pd.DataFrame(all_torsions)
    pairs_df.to_csv(args.output_dir / "pairs.csv", index=False)
    hbonds_df.to_csv(args.output_dir / "hbonds.csv", index=False)
    torsions_df.to_csv(args.output_dir / "torsions.csv", index=False)

    summary = {
        "n_pdbs_attempted": len(pdbs),
        "n_pdbs_succeeded": len(pdbs) - len(failed),
        "n_pdbs_failed": len(failed),
        "n_pairs": len(all_pairs),
        "n_hbonds": len(all_hbonds),
        "n_torsion_observations": len(all_torsions),
        "elapsed_seconds": round(time.time() - t0, 1),
        "failed_pdbs": [{"pdb_id": p, "error": e} for p, e in failed],
    }
    with open(args.output_dir / "extraction_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    print(f"\nWrote feature tables to {args.output_dir}/")
    print(f"  pairs.csv     : {len(all_pairs):,} rows")
    print(f"  hbonds.csv    : {len(all_hbonds):,} rows")
    print(f"  torsions.csv  : {len(all_torsions):,} rows")
    if failed:
        print(f"  Failed PDBs:        {len(failed)} (see extraction_summary.json)")


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
#
# Requires PARSE_CODE_REPO (pair_finder) and PARSE_JSON_DIR (raw per-structure
# JSON, on the data drive). See scripts/data_gen/README.md.
#
# 1. Smoke test on the first 5 PDBs (verifies the pipeline works end-to-end):
#
#       python scripts/data_gen/extract_pair_features.py --limit 5
#
# 2. Full extraction over every JSON-ready PDB in reference_set.json:
#
#       python scripts/data_gen/extract_pair_features.py
#
# 3. Use a different reference set (e.g., a stricter A/B variant):
#
#       python scripts/data_gen/extract_pair_features.py \
#           --reference-set data/reference/reference_set_strict.json \
#           --output-dir   data/reference/features_strict
#
# 4. Point at the raw JSON explicitly (overrides PARSE_JSON_DIR):
#
#       python scripts/data_gen/extract_pair_features.py --json-dir "$PARSE_JSON_DIR"
#
# Output:
#   data/reference/high_quality_features/
#     pairs.csv            one row per base-base pair
#     hbonds.csv           one row per H-bond in a base-base pair
#     torsions.csv         one row per residue × pair it appears in
#     extraction_summary.json  per-run counts + list of failed PDBs
#
# Step 2 (build_distributions.py) consumes these tables to produce the
# per-(bp_type, edge_type, parameter) Z' / ProSco distributions used by
# the scorer at runtime.