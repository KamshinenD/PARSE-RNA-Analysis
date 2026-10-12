"""Time-reversal (RTM) control for the RE-SOLVED cohort.

The PDB-REDO control (time_reversal_controls.py) can match pairs on res_id alone,
because both models are the same entry re-refined. Here the two models are
DIFFERENT PDB entries of the same molecule, so chains are renamed and only the
chain pair recorded in resolve_strict.csv (cw in the worse entry, cb in the better
one) is comparable. Pairs are matched within that chain on residue number AND base
identity, in both orientations, so a numbering offset cannot silently align two
different residues.

Writes data/reference/resolved_again/resolve_pairs_matched.csv (the per-pair table
this cohort never had) and data/validation/time_reversal_resolved.json.
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _rtm import N_BOOT, paired_control, permutation_null  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
STRICT = ROOT / "data/reference/resolved_again/resolve_strict.csv"
sys.path.insert(0, str(ROOT / "scripts" / "data_gen"))
import _engine  # noqa: E402

ENGINE = _engine.find_engine()
CUT = 80.0


def scored_pairs(pdb: str) -> list[dict]:
    """`parse <pdb> --details` -> the pairs that carry a score."""
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as fh:
        out = fh.name
    r = subprocess.run([ENGINE, *_engine.structure_args(pdb), "--details", "--out", out],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"parse failed for {pdb}: {r.stderr.strip()[:160]}")
    d = json.loads(Path(out).read_text())
    Path(out).unlink(missing_ok=True)
    rows = []
    for p in d.get("pairs") or []:
        if p.get("score") is None:
            continue
        rows.append({"res_id1": p["res_id1"], "res_id2": p["res_id2"],
                     "lw_class": (p.get("classification") or {}).get("lw_class"),
                     "score": float(p["score"])})
    return rows


def _split(res_id: str):
    """'A3-C-51' -> (chain, base, number); chain ids can contain '-'-free digits."""
    chain, base, num = res_id.rsplit("-", 2)
    return chain, base, num


def index_by_chain(rows, chain):
    idx = {}
    for p in rows:
        c1, b1, n1 = _split(p["res_id1"])
        c2, b2, n2 = _split(p["res_id2"])
        if c1 != chain or c2 != chain:
            continue
        k = (n1, b1, n2, b2, p["lw_class"])
        idx[k] = p["score"]
        idx[(n2, b2, n1, b1, p["lw_class"])] = p["score"]
    return idx


MATCHED = ROOT / "data/reference/resolved_again/resolve_pairs_matched.csv"


def build_matched() -> pd.DataFrame:
    S = pd.read_csv(STRICT)
    cache: dict[str, list[dict]] = {}
    recs, failed = [], []
    for i, r in S.iterrows():
        try:
            for pdb in (r.worst, r.best):
                if pdb not in cache:
                    cache[pdb] = scored_pairs(pdb)
        except Exception as e:                      # noqa: BLE001
            failed.append((r.rep, r.worst, r.best, str(e)[:80]))
            continue
        iw = index_by_chain(cache[r.worst], str(r.cw))
        ib = index_by_chain(cache[r.best], str(r.cb))
        for k, sw in iw.items():
            if k in ib:
                recs.append({"rep": r.rep, "worst": r.worst, "best": r.best,
                             "num1": k[0], "base1": k[1], "num2": k[2], "base2": k[3],
                             "lw_class": k[4], "score_w": sw, "score_b": ib[k]})
        print(f"[{i + 1:3d}/{len(S)}] {r.rep} {r.worst}->{r.best} "
              f"chains {r.cw}/{r.cb}: {sum(1 for k in iw if k in ib) // 2} matched",
              flush=True)
    df = pd.DataFrame(recs).drop_duplicates(
        subset=["rep", "worst", "best", "num1", "base1", "num2", "base2", "lw_class"])
    # each pair was indexed in both orientations; keep one row per physical pair
    df["key"] = [tuple(sorted([(a, b), (c, d)])) + (lw,) for a, b, c, d, lw
                 in zip(df.num1, df.base1, df.num2, df.base2, df.lw_class)]
    df = df.drop_duplicates(subset=["rep", "worst", "best", "key"]).drop(columns="key")
    df.to_csv(MATCHED, index=False)
    print(f"\nmatched {len(df):,} pairs across {df.rep.nunique()} molecules -> "
          f"{MATCHED.relative_to(ROOT)}")
    if failed:
        print("failed:", failed[:5], f"({len(failed)} total)")
    return df


def main(rebuild: bool = False):
    df = build_matched() if (rebuild or not MATCHED.exists()) else pd.read_csv(MATCHED)

    ctl = paired_control(df, "rep", "score_w", "score_b", CUT)
    ctl["cohort"] = "re-solved (independent re-determination of the same molecule)"
    ctl["source"] = "data/reference/resolved_again/resolve_strict.csv + parse --details"
    ctl["unit"] = "molecule"
    ctl["matching"] = ("chain cw/cb, key = (resnum, base) x2 + lw_class, both orientations")
    ctl["threshold"] = CUT
    ctl["n_pairs_total"] = int(len(df))
    ctl["n_units_total"] = int(df.rep.nunique())
    ctl["estimator"] = ("per-molecule means, paired on molecules present in both "
                        f"directions; bootstrap over molecules, {N_BOOT} resamples")
    ctl["pooled_pair_weighted"] = {
        "forward": float((df[df.score_w < CUT].score_b - df[df.score_w < CUT].score_w).mean()),
        "reverse": float((df[df.score_b < CUT].score_w - df[df.score_b < CUT].score_b).mean()),
        "note": "kept for comparison only: a few ribosomes carry a quarter of the pairs",
    }
    ctl["n_flagged_in_both_models"] = int(((df.score_w < CUT) & (df.score_b < CUT)).sum())
    ctl["permutation_null"] = permutation_null(df, "rep", "score_w", "score_b", CUT)
    # sensitivity: molecules whose two records barely agree behave like that null
    ctl["sensitivity_r_floor_0.2"] = paired_control(df, "rep", "score_w", "score_b",
                                                    CUT, r_floor=0.2)
    (ROOT / "data/validation/time_reversal_resolved.json").write_text(
        json.dumps(ctl, indent=1) + "\n")
    a = ctl["asymmetry"]; t = ctl["implied_true_improvement"]
    print(f"\nmolecules paired: {ctl['n_units']}  (forward-only {ctl['n_units_forward_only']}, "
          f"reverse-only {ctl['n_units_reverse_only']})")
    for d in ("forward", "reverse"):
        b = ctl[d]
        print(f"  {d:8s} {b['mean']:+6.2f}  CI [{b['ci95'][0]:+.2f}, {b['ci95'][1]:+.2f}]")
    print(f"  asymmetry {a['mean']:+6.2f}  CI [{a['ci95'][0]:+.2f}, {a['ci95'][1]:+.2f}]  "
          f"P(<=0)={a['p_le_0']:.3f}  ({a['units_with_forward_greater']}/{ctl['n_units']} molecules)")
    print(f"  implied true improvement {t['mean']:+.2f}  "
          f"CI [{t['ci95'][0]:+.2f}, {t['ci95'][1]:+.2f}]")
    s2 = ctl["sensitivity_r_floor_0.2"]
    print(f"  sensitivity (drop molecules with r < 0.2): {s2['n_units']} molecules, "
          f"asymmetry {s2['asymmetry']['mean']:+.2f} "
          f"[{s2['asymmetry']['ci95'][0]:+.2f}, {s2['asymmetry']['ci95'][1]:+.2f}]")
    print(f"  permutation null: forward {ctl['permutation_null']['forward_mean']:+.2f}, "
          f"reverse {ctl['permutation_null']['reverse_mean']:+.2f} (symmetric by construction)")


if __name__ == "__main__":
    sys.exit(main(rebuild="--rebuild" in sys.argv))
