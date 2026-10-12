"""Referee votes on FastRNA-DSSR class conflicts, restricted to the conflicts that
Supplementary Figure S2A also counts.

Every class conflict is filed twice in the differences inventory
(data/comparison/{class}.json): once under DSSR's class, as
``python_assigns:<FastRNA class>``, and once under FastRNA's class, as
``dssr_assigns:<DSSR class>``. The published referee tally in
validation/annotator_consensus.json votes on the first form (15,356 records);
S2A counts the second form with DSSR's incomplete classes (a '.' edge) dropped
(11,756). The two overlap in 10,194 pairs: the rest are mostly FastRNA's ambiguous
X|Y calls, whose label is not one of the 12 class names and so is never filed
under a FastRNA class. This script re-runs the same vote on the shared 10,194, so
S2A and S2B describe the same conflicts.

The vote is the one in adjudicate_consensus.py (class_disagreement branch), on
_vendor/annotators.py: direction-agnostic; a referee class equal to the disputed
DSSR class or its flip counts for DSSR, one equal to a remaining FastRNA candidate
counts for FastRNA; records whose every FastRNA candidate is direction-equivalent
to the disputed class are degenerate and excluded. Before writing anything the
script re-runs the full population and refuses to continue unless it reproduces
annotator_consensus.json exactly.

Writes data/validation/annotator_consensus_shared.json. The original
annotator_consensus.json is not touched (Supplementary Figure S2D also reads it).
"""
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _vendor import annotators  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
REFEREES = ("fr3d", "rnaview", "barnaba", "rnapolis", "bpnet")
VOTING_TOOLS = ("rnaview", "barnaba", "rnapolis", "bpnet")
OUT = ROOT / "data" / "validation" / "annotator_consensus_shared.json"


def s1a_conflicts():
    """Pairs S2A counts as a class conflict: filed under FastRNA's class with a
    complete DSSR class (no '.' edge)."""
    keys = set()
    for lw in annotators.loaders.LW_COMPARISON_CLASSES:
        rows = json.loads(annotators.paths.data_file(f"comparison/{lw}.json").read_text())
        for r in rows["differences"]:
            reason = r.get("reason") or ""
            if (reason.startswith("dssr_assigns:")
                    and "." not in reason.split(":", 1)[1]
                    and "." not in (r.get("dssr_lw_class") or "")):
                keys.add((r["pdb_id"], frozenset((r["res_id1"], r["res_id2"]))))
    return keys


def class_records(keep=None):
    """Deduplicated class-disagreement records, as annotators.disagreement_records
    yields them, optionally restricted to the pair keys in ``keep``."""
    seen = set()
    for lw in annotators.loaders.LW_COMPARISON_CLASSES:
        rows = json.loads(annotators.paths.data_file(f"comparison/{lw}.json").read_text())
        for row in rows["differences"]:
            key = (row["pdb_id"], frozenset((row["res_id1"], row["res_id2"])), row["source"])
            if key in seen:
                continue
            seen.add(key)
            rec = annotators._record(row, lw)
            if rec.bucket != "class_disagreement":
                continue
            if keep is not None and key[:2] not in keep:
                continue
            yield rec


def tally(records, panels_by_tool):
    per_tool, degenerate, n = Counter(), 0, 0
    for r in records:
        found = {t: {c for c in annotators.lookup(p[r.pdb_id], r.res1, r.res2)
                     if not c.startswith("n")}
                 for t, p in panels_by_tool.items() if r.pdb_id in p}
        if not found:
            continue
        sides = annotators.vote_sides(r.disputed_class, r.parse_candidates)
        if sides is None:
            degenerate += 1
            continue
        n += 1
        dssr_side, parse_side = sides
        for tool, classes in found.items():
            per_tool[(tool, "covered")] += 1
            if classes & parse_side:
                per_tool[(tool, "votes_parse")] += 1
            if classes & dssr_side:
                per_tool[(tool, "votes_dssr")] += 1
    return {t: {k: per_tool[(t, k)] for k in ("covered", "votes_parse", "votes_dssr")}
            for t in REFEREES}, degenerate, n


def main():
    panels = {"fr3d": annotators.fr3d_pairs(),
              **{t: annotators.tool_pairs(t) for t in VOTING_TOOLS}}

    full, _, _ = tally(class_records(), panels)
    published = json.loads((ROOT / "data/validation/annotator_consensus.json").read_text())
    pub = published["per_annotator"]
    for t in REFEREES:
        for k in ("covered", "votes_parse", "votes_dssr"):
            want = pub.get(f"class_disagreement|{t}|{k}", 0)
            if full[t][k] != want:
                sys.exit(f"does not reproduce annotator_consensus.json: {t} {k} "
                         f"{full[t][k]} != {want}")
    print("full population reproduces annotator_consensus.json exactly")

    keep = s1a_conflicts()
    shared, degenerate, n = tally(class_records(keep), panels)
    out = {
        "population": "class conflicts counted in S2A: filed under FastRNA's class with "
                      "a complete DSSR class, and voted from DSSR's side as in "
                      "annotator_consensus.json",
        "n_s1a_conflicts": len(keep),
        "n_records_voted": n,
        "degenerate_excluded": degenerate,
        "referees": list(REFEREES),
        "per_annotator": {f"class_disagreement|{t}|{k}": v
                          for t, d in shared.items() for k, v in d.items()},
        "provenance": "same July differences inventory and referee files as "
                      "annotator_consensus.json; vote logic from "
                      "adjudicate_consensus.py (see module docstring)",
    }
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    for t in REFEREES:
        d = shared[t]
        print(f"  {t:9s} covered {d['covered']:6,d}  FastRNA "
              f"{100 * d['votes_parse'] / d['covered']:5.1f}%  DSSR "
              f"{100 * d['votes_dssr'] / d['covered']:5.1f}%")
    print("wrote", OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
