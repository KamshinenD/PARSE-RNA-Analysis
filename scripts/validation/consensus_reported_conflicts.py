"""Referee consensus on FastRNA-DSSR class conflicts, restricted to pairs FastRNA
actually reports.

The comparison inventory (data/comparison/{class}.json) files a DSSR pair as a class
conflict ("python_assigns:<class>") whenever FastRNA classified a candidate for it
differently -- including candidates FastRNA then did NOT select, so never reported.
3,817 of the class-conflict records in the referee-panel structures are such
unreported candidates. Those are not FastRNA's calls, so this re-runs the same
adjudication with them removed: a class conflict is kept only if FastRNA reports the
pair (it is present in comparison/parse_pairs.json, the current engine's output).
Every other bucket (rejections, FastRNA-only pairs) is left exactly as it was.

Ambiguous FastRNA calls (X|Y) are scored by what FastRNA actually reported: a
referee that reports ANY of FastRNA's candidates supports FastRNA. The published
rule dropped a candidate whenever it equalled DSSR's class (or its flip), so a
referee agreeing with that shared candidate counted for DSSR only. Here it counts
for both tools, since both named that class; such records become ties rather than
DSSR wins. Single-class calls are scored exactly as before.

The adjudication is adjudicate_consensus.py unchanged; before writing anything this
script re-runs it on the full inventory and refuses to continue unless it reproduces annotator_consensus.json.

Writes data/validation/annotator_consensus_reported.json; annotator_consensus.json is
not touched.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adjudicate_consensus as ac  # noqa: E402
from _vendor import annotators  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PUBLISHED = ROOT / "data" / "validation" / "annotator_consensus.json"
OUT = ROOT / "data" / "validation" / "annotator_consensus_reported.json"


def summary(tallies, tools, dest):
    ac.OUT = dest
    ac.write_summary(tallies, tools)
    return json.loads(dest.read_text())


def main():
    fr3d = annotators.fr3d_pairs()
    tools = ac.load_available(list(ac.VOTING_TOOLS))
    scratch = OUT.with_suffix(".check.json")

    full = summary(ac.run(tools, fr3d, None), tools, scratch)
    scratch.unlink()
    published = json.loads(PUBLISHED.read_text())
    if {k: v for k, v in full.items() if k != "fr3d_local_vs_published"} != \
            {k: v for k, v in published.items() if k != "fr3d_local_vs_published"}:
        sys.exit("full run does not reproduce annotator_consensus.json; not writing")
    print("full run reproduces annotator_consensus.json")

    reported = annotators.parse_pairs()
    original = annotators.disagreement_records
    dropped = {"n": 0}

    def reported_only():
        for r in original():
            # only structures with FastRNA output in parse_pairs.json (the referee
            # panel) can be checked; conflicts elsewhere have no referee coverage
            # and pass through untouched rather than being counted as unreported
            if (r.bucket == "class_disagreement" and r.pdb_id in reported
                    and not annotators.lookup(reported[r.pdb_id], r.res1, r.res2)):
                dropped["n"] += 1
                continue
            yield r

    published_sides = annotators.vote_sides

    def vote_sides(disputed, candidates):
        # ambiguous call: every FastRNA candidate (and its flip) supports FastRNA
        if len(candidates) > 1:
            dssr = {disputed, annotators.flip(disputed)}
            return dssr, {v for c in candidates for v in (c, annotators.flip(c))}
        return published_sides(disputed, candidates)

    annotators.disagreement_records = reported_only
    annotators.vote_sides = vote_sides
    try:
        out = summary(ac.run(tools, fr3d, None), tools, OUT)
    finally:
        annotators.disagreement_records = original
        annotators.vote_sides = published_sides

    out["provenance_note"] = (
        "same comparison inventory and referee files as annotator_consensus.json (the "
        "inventory was re-run with the current C++ engine on 2026-10-01 and is "
        "identical row for row); class conflicts restricted to pairs FastRNA reports "
        f"(present in comparison/parse_pairs.json): {dropped['n']:,} unreported-candidate "
        "records removed; ambiguous FastRNA calls credited when a referee reports any "
        "of FastRNA's candidates (a candidate shared with DSSR credits both); other "
        "buckets unchanged")
    out["class_disagreement_unreported_removed"] = dropped["n"]
    OUT.write_text(json.dumps(out, indent=1))
    pops = out["class_disagreement_populations"]
    for p in ("clean", "shared_candidate", "pooled"):
        print(f"  {p:17s} {pops[p]['verdicts']}")
    print(f"removed {dropped['n']:,} unreported class-conflict records; wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
