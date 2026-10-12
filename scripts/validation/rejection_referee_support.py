"""How many of the five referees report the pairs FastRNA rejects, split by reason.

annotator_consensus.json tallies referee support for every rejection pooled
(``parse_rejections_n_annotators_reporting``, 80,116 records). Supplementary Figure
S2A shows only one rejection reason -- the pair failed candidate generation
(``not_candidate``) -- and leaves the selection-stage reasons (no strong H-bond, not
selected, stacking) out, so S2's referee panel needs the same split. This script
tallies each reason separately, with the referee panel and lookup that
adjudicate_consensus.py uses, and refuses to write unless the reasons summed back
together reproduce the published pooled tally exactly.

Writes data/validation/rejection_referee_support.json.
"""
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adjudicate_consensus as ac  # noqa: E402
from _vendor import annotators  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "validation" / "rejection_referee_support.json"


def main():
    panels = {"fr3d": annotators.fr3d_pairs(),
              **ac.load_available(list(ac.VOTING_TOOLS))}
    by_reason = defaultdict(Counter)          # reason -> (panel_size, n_reporting)
    for r in annotators.disagreement_records():
        if r.bucket != "rejected_by_parse":
            continue
        found = {t: ac.real_classes(p[r.pdb_id], r) for t, p in panels.items()
                 if r.pdb_id in p}
        if not found:
            continue
        n_rep = sum(1 for classes in found.values() if classes)
        by_reason[r.reason][(len(found), n_rep)] += 1

    pooled = sum(by_reason.values(), Counter())
    published = json.loads((ROOT / "data/validation/annotator_consensus.json").read_text())
    pub = published["parse_rejections_n_annotators_reporting"]
    mine = {f"panel_{p}": {} for p, _ in pooled}
    for (p, n), v in pooled.items():
        mine[f"panel_{p}"][str(n)] = v
    norm = lambda d: {k: {str(a): b for a, b in v.items()} for k, v in d.items()}
    if norm(mine) != norm(pub):
        sys.exit("pooled tally does not reproduce annotator_consensus.json; not writing")
    print("pooled tally reproduces annotator_consensus.json")

    out = {"note": "referees (of the five) reporting each FastRNA-rejected pair, by "
                   "rejection reason; panel_5 = records all five referees cover",
           "by_reason": {}}
    for reason, c in sorted(by_reason.items()):
        block = defaultdict(dict)
        for (p, n), v in sorted(c.items()):
            block[f"panel_{p}"][str(n)] = v
        out["by_reason"][reason] = block
        n5 = sum(block.get("panel_5", {}).values())
        print(f"  {reason:16s} panel_5 n={n5:6,d}  " + " ".join(
            f"{k}:{100 * v / n5:4.1f}%" for k, v in sorted(block.get("panel_5", {}).items())))
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print("wrote", OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
