"""Score original (RCSB) vs PDB-REDO re-refined coordinates with the C++ engine.

For each structure that has a cached PDB-REDO file (run fetch_pdb_redo.py first),
run the production ``parse`` binary twice:

    * original  ->  parse <PDBID>                    (parse fetches from RCSB)
    * PDB-REDO  ->  parse redo_files/<PDB>/<PDB>.pdb  (fetched _final.pdb)

and write both validation files (raw scores only; deltas computed in-notebook):

    data/validation/pdb_redo_comparison_all.json    (structure-level aggregates)
    data/validation/pdb_redo_pairs_comparison.json  (per matched base pair)

Build the engine and set PARSE_ENGINE (or put `parse` on PATH) — see the repo
README. The validation set is whatever REDO files are present under --redo-cache.
"""

import argparse
import datetime
import json
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts" / "data_gen"))
import _engine  # noqa: E402

REDO_CACHE = REPO_ROOT / "data" / "validation" / "redo_files"
OUT_ALL = REPO_ROOT / "data" / "validation" / "pdb_redo_comparison_all.json"
OUT_PAIRS = REPO_ROOT / "data" / "validation" / "pdb_redo_pairs_comparison.json"

_ENGINE: str | None = None
_REDO_CACHE: Path = REDO_CACHE


def _init_worker(engine: str, redo_cache: Path) -> None:
    global _ENGINE, _REDO_CACHE
    _ENGINE = engine
    _REDO_CACHE = redo_cache


def run_parse(arg: str) -> dict | None:
    """Run `parse` on a PDB id or coordinate file; return its JSON (stdout)."""
    try:
        out = subprocess.run([_ENGINE, arg], capture_output=True, text=True, timeout=900)
    except Exception:  # noqa: BLE001
        return None
    if out.returncode != 0 or not out.stdout.strip():
        return None
    try:
        return json.loads(out.stdout)
    except json.JSONDecodeError:
        return None


def extract(d: dict) -> tuple[dict, list[dict]]:
    """(structure-level scores, per-pair records) from one parse output.

    Only unambiguously classified pairs carry a numeric score; ambiguous
    ("X|Y") pairs have score None and are counted as skipped.
    """
    pairs = d.get("pairs", [])
    scored = [p for p in pairs if isinstance(p.get("score"), (int, float))]
    issues: Counter = Counter()
    for p in scored:
        for iss in (p.get("issues") or []):
            issues[iss] += 1
    puck = sorted({u.get("res_id") for p in pairs
                   for u in (p.get("unusual_pucker") or []) if u.get("res_id")})
    struct = {
        "pairs_score": d.get("pairs_score"),
        "backbone_suiteness": d.get("backbone_suiteness"),
        "n_pairs": d.get("n_pairs"),
        "n_residues_scored": len(d.get("backbone_residues") or []),
        "skipped_pairs": len(pairs) - len(scored),
        "issue_summary": dict(issues),
        "tier_summary": d.get("tier_summary") or {},
        "n_unusual_pucker": len(puck),
        "unusual_pucker_res": puck,
    }
    pair_recs = [{
        "res_id1": p["res_id1"],
        "res_id2": p["res_id2"],
        "lw_class": (p.get("classification") or {}).get("lw_class"),
        "bp_type": p.get("sequence"),
        "score": p["score"],
        "issues": p.get("issue_details") or [],
    } for p in scored]
    return struct, pair_recs


def match_pairs(orig: list[dict], redo: list[dict]) -> list[dict]:
    """Match by (res_id1, res_id2, lw_class), trying both orientations."""
    idx: dict[tuple, dict] = {}
    for p in redo:
        idx[(p["res_id1"], p["res_id2"], p["lw_class"])] = p
        idx[(p["res_id2"], p["res_id1"], p["lw_class"])] = p
    out = []
    for op in orig:
        rp = idx.get((op["res_id1"], op["res_id2"], op["lw_class"]))
        if rp is None:
            continue
        out.append({
            "res_id1": op["res_id1"], "res_id2": op["res_id2"],
            "lw_class": op["lw_class"], "bp_type": op["bp_type"],
            "orig_score": op["score"], "redo_score": rp["score"],
            "orig_issues": op["issues"], "redo_issues": rp["issues"],
        })
    return out


def process_one(pdb: str) -> dict:
    """Score original + REDO for one PDB; return the combined record."""
    do = run_parse(pdb)
    if do is None:
        return {"status": "fail", "failure": {"pdb": pdb, "stage": "original", "error": "parse failed"}}
    redo_file = _REDO_CACHE / pdb / f"{pdb}.pdb"
    if not redo_file.exists():
        return {"status": "fail", "failure": {"pdb": pdb, "stage": "redo", "error": "no cached redo pdb"}}
    dr = run_parse(str(redo_file))
    if dr is None:
        return {"status": "fail", "failure": {"pdb": pdb, "stage": "redo", "error": "parse failed"}}
    o_struct, o_pairs = extract(do)
    r_struct, r_pairs = extract(dr)
    matched = match_pairs(o_pairs, r_pairs)
    for m in matched:
        m["pdb_id"] = pdb
    return {"status": "ok",
            "result": {"pdb_id": pdb, "original": o_struct, "redo": r_struct},
            "matched": matched,
            "n_unmatched": len(o_pairs) - len(matched)}


def discover_cached(redo_cache: Path) -> list[str]:
    """PDB IDs with a cached REDO coordinate file (<PDB>/<PDB>.pdb)."""
    out = []
    for d in sorted(redo_cache.iterdir()):
        if d.is_dir() and (d / f"{d.name}.pdb").exists():
            out.append(d.name.upper())
    return out


def _stat(xs):
    xs_s = sorted(xs)
    return {"mean": round(sum(xs) / len(xs), 3), "median": round(xs_s[len(xs_s) // 2], 3),
            "min": round(min(xs), 3), "max": round(max(xs), 3)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--redo-cache", type=Path, default=REDO_CACHE)
    ap.add_argument("--pdb-list", type=Path, default=None,
                    help="explicit id list (default: every cached REDO structure)")
    ap.add_argument("--engine", default=None, help="path to `parse` (else PARSE_ENGINE / PATH)")
    ap.add_argument("--out-all", type=Path, default=OUT_ALL)
    ap.add_argument("--out-pairs", type=Path, default=OUT_PAIRS)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    engine = _engine.find_engine(args.engine)
    if args.pdb_list:
        pdbs = [p.strip().upper() for p in args.pdb_list.read_text().splitlines()
                if p.strip() and not p.startswith("#")]
    else:
        pdbs = discover_cached(args.redo_cache)
    if args.limit:
        pdbs = pdbs[:args.limit]
    print(f"Engine: {engine}\nRescoring {len(pdbs)} cached structures (workers={args.workers})\n", flush=True)

    results, all_pairs, failures = [], [], []
    n_unmatched = 0
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers,
                             initializer=_init_worker, initargs=(engine, args.redo_cache)) as ex:
        for i, rec in enumerate(ex.map(process_one, pdbs), 1):
            if rec["status"] == "ok":
                results.append(rec["result"])
                all_pairs.extend(rec["matched"])
                n_unmatched += rec["n_unmatched"]
            else:
                failures.append(rec["failure"])
            if i % 25 == 0 or i == len(pdbs):
                el = time.time() - t0
                rate = i / el if el else 0
                print(f"  [{i}/{len(pdbs)}] paired={len(results)} pairs={len(all_pairs):,} "
                      f"failed={len(failures)} | {rate:.1f}/s", flush=True)
    elapsed = round(time.time() - t0, 1)

    ps = {"original": _stat([r["original"]["pairs_score"] for r in results if r["original"]["pairs_score"] is not None]),
          "redo": _stat([r["redo"]["pairs_score"] for r in results if r["redo"]["pairs_score"] is not None])} if results else {}
    all_payload = {
        "schema_version": "3.0",
        "generated_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "method": {
            "_documentation": [
                "Structure-level scores on original (RCSB) vs PDB-REDO coordinates,",
                "both scored by the production C++ `parse` binary (find+classify+score).",
                "Schema 3.0: no `overall`; backbone reported as `backbone_suiteness` (0-1).",
                "Raw scores only; deltas computed in-notebook.",
            ],
            "scorer": "cpp:parse",
        },
        "summary": {"n_attempted": len(pdbs), "n_paired": len(results),
                    "n_failed": len(failures), "pairs_score": ps},
        "results": results, "failures": failures,
    }
    args.out_all.write_text(json.dumps(all_payload, indent=2, ensure_ascii=False))

    if all_pairs:
        o = [p["orig_score"] for p in all_pairs]
        r = [p["redo_score"] for p in all_pairs]
        deltas = [b - a for a, b in zip(o, r)]
        wins = sum(1 for d in deltas if d > 0)
        pair_summary = {
            "n_structures": len(results), "n_pairs_matched": len(all_pairs),
            "n_unmatched": n_unmatched, "n_failed": len(failures),
            "elapsed_seconds": elapsed,
            "orig_score": _stat(o), "redo_score": _stat(r), "delta": _stat(deltas),
            "redo_wins": wins, "redo_wins_pct": round(100 * wins / len(all_pairs), 2),
        }
    else:
        pair_summary = {"n_structures": 0, "n_pairs_matched": 0}
    pairs_payload = {
        "schema_version": "2.0",
        "generated_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "method": {
            "_documentation": [
                "Per matched base pair: score on original vs PDB-REDO coordinates,",
                "both from the production C++ `parse` binary.",
                "Matched by (res_id1, res_id2, lw_class); both orientations tried.",
            ],
            "scorer": "cpp:parse",
        },
        "summary": pair_summary, "pairs": all_pairs, "failures": failures,
    }
    args.out_pairs.write_text(json.dumps(pairs_payload, indent=2, ensure_ascii=False))

    print(f"\nWrote {args.out_all.name}: {len(results)} paired, {len(failures)} failed")
    print(f"Wrote {args.out_pairs.name}: {len(all_pairs):,} matched pairs ({n_unmatched:,} unmatched)")
    if results:
        print(f"  pairs_score  original mean={ps['original']['mean']:.2f}  redo mean={ps['redo']['mean']:.2f}")
    if all_pairs:
        print(f"  per-pair Δ mean={pair_summary['delta']['mean']:+.2f}  "
              f"redo wins {pair_summary['redo_wins']:,}/{len(all_pairs):,} "
              f"({pair_summary['redo_wins_pct']:.1f}%)")


if __name__ == "__main__":
    main()
