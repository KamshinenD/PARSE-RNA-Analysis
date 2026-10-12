"""Referee coverage of FastRNA's rejections and its unique calls, current engine.

The published tally in validation/annotator_consensus.json was computed from the
July differences inventory (data/comparison/{class}.json, 2026-07-22); its own
provenance note records that 24.9% of those disagreement pairs are absent from the
current engine output and that 1,484 rejection records are now called. This script
recomputes the two coverage panels directly from the CURRENT pair files, so the
figure describes the engine that ships.

Two populations, both restricted to structures every tool processed successfully:

  rejected  - DSSR reports the pair, FastRNA does not
  unique    - FastRNA reports the pair, DSSR does not

For each pair we count how many of the five independent referees (published FR3D,
RNAView, barnaba, RNApolis, BPNET) also report it. Matching is direction-agnostic
on (chain, residue number) and ignores the class label: the question here is
whether the interaction is seen at all, not how it is classified.
"""
import csv
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "data" / "comparison"
REFEREES = ["fr3d", "rnaview", "barnaba", "rnapolis", "bpnet"]


def load_tool(name):
    """{pdb: set of direction-agnostic pair keys}, only for structures it handled."""
    if name == "fr3d":                       # published FR3D: one csv per structure
        out = {}
        for f in sorted((COMP / "fr3d").glob("*.csv")):
            keys = set()
            for row in csv.reader(f.open()):
                if len(row) < 3:
                    continue
                a, b = row[0].split("|"), row[2].split("|")
                if len(a) < 5 or len(b) < 5:
                    continue
                keys.add(frozenset({(a[2], a[4]), (b[2], b[4])}))
            out[f.stem.upper()] = keys
        return out
    raw = json.loads((COMP / f"{name}_pairs.json").read_text())
    return {pdb.upper(): {frozenset({(p[0], str(p[1])), (p[2], str(p[3]))})
                          for p in (e.get("pairs") or [])}
            for pdb, e in raw.items() if e.get("status") == "ok"}


def main():
    tools = {n: load_tool(n) for n in ["parse", "dssr"] + REFEREES}
    for n, d in tools.items():
        print(f"  {n:9s} {len(d):5,d} structures", flush=True)

    # a structure counts only where every tool succeeded, so the panel size is 5
    shared = set.intersection(*(set(d) for d in tools.values()))
    print(f"\n{len(shared):,} structures processed by all {len(tools)} tools")

    out = {"source": "data/comparison/*.json + fr3d/ (current engine run)",
           "referees": REFEREES, "n_structures": len(shared),
           "matching": "direction-agnostic on (chain, resnum); class label ignored",
           "populations": {}}
    for label, a, b in (("rejected_by_fastrna", "dssr", "parse"),
                        ("fastrna_only", "parse", "dssr")):
        hist, total = Counter(), 0
        for pdb in shared:
            only = tools[a][pdb] - tools[b][pdb]
            if not only:
                continue
            refs = [tools[r][pdb] for r in REFEREES]
            for key in only:
                hist[sum(key in r for r in refs)] += 1
                total += 1
        counts = {str(k): hist.get(k, 0) for k in range(6)}
        pct = {k: 100 * v / total for k, v in counts.items()}
        out["populations"][label] = {
            "n": total, "counts": counts,
            "pct_by_n_referees": {k: round(v, 2) for k, v in pct.items()},
            "pct_ge_3": round(sum(pct[str(k)] for k in (3, 4, 5)), 2),
            "pct_le_1": round(sum(pct[str(k)] for k in (0, 1)), 2)}
        print(f"\n{label}: {total:,} pairs")
        print("   referees reporting: " +
              "  ".join(f"{k}:{pct[str(k)]:.0f}%" for k in range(6)))
        print(f"   >= 3 referees: {out['populations'][label]['pct_ge_3']:.0f}%   "
              f"<= 1 referee: {out['populations'][label]['pct_le_1']:.0f}%")
    dest = ROOT / "data/validation/referee_coverage_current.json"
    dest.write_text(json.dumps(out, indent=1) + "\n")
    print("\nwrote", dest.relative_to(ROOT))


if __name__ == "__main__":
    sys.exit(main())
