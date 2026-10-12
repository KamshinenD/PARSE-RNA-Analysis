"""Build the re-solved cohort: the same RNA deposited twice at different resolution.

    python scripts/validation/build_resolved_cohort.py            # -> resolve_strict.csv
    python scripts/validation/build_resolved_cohort.py --out x.csv

Writes data/reference/resolved_again/resolve_strict.csv (one row per molecule), which
time_reversal_resolved.py and the Figure 5 / Supplementary S7 tables are built from.

Selection (Methods, "Re-determination from better data"):
  1. Members come from the BGSU equivalence classes: the `every` column of
     uniqueRNAS.csv lists every deposited chain of each class. A member is kept if
     it has a known resolution and engine output for that chain.
  2. A class is used only if it has >= 2 such members, all from one experimental
     method (no X-ray vs cryo-EM comparisons).
  3. Every pair of members whose resolutions differ by >= MIN_GAP is a candidate;
     the lower-resolution entry is the original ("worst"), the other the
     re-determined model ("best").
  4. Candidates from all classes are taken greedily in order of decreasing gap,
     at most one per class, and no PDB entry is ever used twice.
  5. A candidate is accepted only if the two chains share >= MIN_SHARED residue
     numbers and >= SEQ_MIN % of those positions carry the same nucleotide
     (A/C/G/U; position-by-position at matching residue numbers, not an alignment).
Pairs are therefore chosen on resolution alone, never on a FastRNA score.

Per molecule it then records the four measures over matched objects:
  base pairs   matched within the chosen chain on (residue 1, residue 2, LW class);
               mean score of those in the Review tier in the original model
  suites       mean suiteness and % unassignable (`!!`) over matched residues
  riboses      count of delta/Pperp pucker outliers over matched residues

Inputs, all under data/: the engine output over the non-redundant set and over the
extra members of each equivalence class (data/reference/engine_archive/, see
_archive.py), uniqueRNAS.csv and the PDB metadata. It reproduces the shipped table
exactly apart from `nt`, which reads the recounted nucleotide totals
(scripts/reference/count_nucleotides.py).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _archive as A  # noqa: E402

OUT = A.DATA / "reference" / "resolved_again" / "resolve_strict.csv"

MIN_GAP = 0.5      # A of resolution between the two depositions
MIN_SHARED = 30    # residue numbers the two chains must share
SEQ_MIN = 95.0     # % identical nucleotides at the shared residue numbers


def split_res(s: pd.Series) -> tuple[pd.Series, pd.Series]:
    """'A-G-12' -> (chain 'A', number 12)."""
    p = s.str.split("-")
    return p.str[0], pd.to_numeric(p.str[-1], errors="coerce")


def residues() -> pd.DataFrame:
    """One row per (entry, chain, residue number), RNA bases only."""
    old = A.residues_orig()[["pdb_id", "chain", "num", "base_type", "suiteness",
                             "is_outlier", "pucker_outlier"]]
    new = A.residues_expand()
    new["chain"], new["num"] = split_res(new.res_id)
    S = pd.concat([old, new[old.columns]], ignore_index=True)
    S["is_outlier"] = S.is_outlier.fillna(False).astype(bool)
    S["pucker_outlier"] = S.pucker_outlier.fillna(False).astype(bool)
    S = S.dropna(subset=["num"]).drop_duplicates(["pdb_id", "chain", "num"])
    return S[S.base_type.isin(["A", "C", "G", "U"])]


def pairs() -> pd.DataFrame:
    """One row per selected base pair, with chain and residue number of each side."""
    old = A.pairs_orig()
    new = A.pairs_expand()
    P = pd.concat([old, new[old.columns]], ignore_index=True)
    P["c1"], P["n1"] = split_res(P.res_id1)
    P["c2"], P["n2"] = split_res(P.res_id2)
    return P.dropna(subset=["n1", "n2"])


def candidates(byc: dict, have: set, res: dict, meth: dict) -> pd.DataFrame:
    """Every same-method (worst, best) pair inside a BGSU class with gap >= MIN_GAP."""
    u = A.uniques()
    rows = []
    for ci, r in enumerate(u.itertuples()):
        members = {}
        for m in str(r.every).split(","):          # "PDB|model|chain"
            if "|" not in m:
                continue
            parts = m.strip().split("|")
            pid, chain = parts[0].upper(), parts[-1]
            if pid in have and (pid, chain) in byc and res.get(pid) is not None \
                    and res.get(pid) == res.get(pid):   # resolution known (not NaN)
                members[pid] = chain
        methods = {meth.get(p) for p in members}
        if len(members) < 2 or len(methods) != 1:
            continue
        method = next(iter(methods))
        for b in members:
            for w in members:
                if res[w] - res[b] >= MIN_GAP:
                    rows.append({"mol": ci, "rep": r.unique_pdb_id, "worst": w, "best": b,
                                 "cw": members[w], "cb": members[b],
                                 "gap": res[w] - res[b], "method": method})
    return pd.DataFrame(rows).sort_values("gap", ascending=False)


def select(cands: pd.DataFrame, byc: dict) -> list:
    """Greedy, widest gap first: one pair per class, no entry reused, same molecule."""
    used_pdb, used_mol, picked = set(), set(), []
    for r in cands.itertuples():
        if r.mol in used_mol or r.worst in used_pdb or r.best in used_pdb:
            continue
        M = byc[(r.worst, r.cw)].merge(byc[(r.best, r.cb)], on="num", suffixes=("_w", "_b"))
        if len(M) < MIN_SHARED:
            continue
        if 100 * (M.base_type_w == M.base_type_b).mean() < SEQ_MIN:
            continue
        used_mol.add(r.mol)
        used_pdb.update((r.worst, r.best))
        picked.append((r, M))
    return picked


def summarise(picked: list, byp: dict, n_nt: dict) -> pd.DataFrame:
    rows = []
    rv = pd.DataFrame()       # Review pairs of the last molecule that had pair tables
    for r, M in picked:
        pa, pb = byp.get((r.worst, r.cw)), byp.get((r.best, r.cb))
        dsc = sw = sb = np.nan
        if pa is not None and pb is not None:
            pa, pb = pa[pa.c2 == r.cw], pb[pb.c2 == r.cb]       # both residues in the chain
            if len(pa) and len(pb):
                mp = pa.merge(pb, on=["n1", "n2", "lw_class"], suffixes=("_w", "_b"))
                rv = mp[mp.tier_w == "Review"].dropna(subset=["score_w", "score_b"])
                if len(rv) >= 1:
                    dsc = (rv.score_b - rv.score_w).mean()
                    sw, sb = rv.score_w.mean(), rv.score_b.mean()
        rows.append({
            "rep": r.rep, "worst": r.worst, "best": r.best, "gap": r.gap,
            "cw": r.cw, "cb": r.cb,
            # As in the shipped table: when one side has no pair table, this carries
            # the previous molecule's count; such rows have no score (score_w is NaN).
            "n_review": 0 if pa is None else len(rv),
            "method": r.method, "n_res": len(M), "nt": n_nt.get(r.best),
            "d_suit": M.suiteness_b.mean() - M.suiteness_w.mean(),
            "d_bang": 100 * (M.is_outlier_b.mean() - M.is_outlier_w.mean()),
            "d_puck": 100 * (M.pucker_outlier_b.mean() - M.pucker_outlier_w.mean()),
            "d_score": dsc, "score_w": sw, "score_b": sb,
            "puck_n_w": int(M.pucker_outlier_w.sum()), "puck_n_b": int(M.pucker_outlier_b.sum()),
            "bang_n_w": int(M.is_outlier_w.sum()), "bang_n_b": int(M.is_outlier_b.sum()),
            "suit_w": M.suiteness_w.mean(), "suit_b": M.suiteness_b.mean(),
            "bang_w": 100 * M.is_outlier_w.mean(), "bang_b": 100 * M.is_outlier_b.mean(),
        })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    S = residues()
    byc = {k: v for k, v in S.groupby(["pdb_id", "chain"])}
    byp = {k: v for k, v in pairs().groupby(["pdb_id", "c1"])}
    md = A.metadata().set_index("pdb_id")
    picked = select(candidates(byc, set(S.pdb_id.unique()), md.resolution.to_dict(),
                               md.method.to_dict()), byc)
    T = summarise(picked, byp, md.n_nt.to_dict())
    T.to_csv(args.out, index=False)

    ids = set(T.worst) | set(T.best)
    print(f"{len(T)} molecules, {len(ids)} distinct entries "
          f"({T.method.value_counts().to_dict()}); gap {T.gap.min():.2f}-{T.gap.max():.2f} A, "
          f"median {T.gap.median():.2f}; shared residues >= {T.n_res.min()} -> {args.out}")


if __name__ == "__main__":
    main()
