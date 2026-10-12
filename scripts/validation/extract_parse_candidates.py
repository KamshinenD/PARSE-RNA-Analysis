"""FastRNA's identified-and-classified candidates, before selection, for every structure.

Figure 2 and S3 use these candidates to credit discarded pairs when the comparison
tool reports the same pair. The export itself retains the full candidate set.
data/comparison/parse_pairs.json holds only the selected output, so this runs the C++
engine's `dump-selection` command, which tags every candidate with the production
selection's own disposition, and keeps the identified ones:

    selected     in the output
    competition  valid and eligible, displaced by a competing pair   (not_selected)
    ineligible   fails class_valid, no strong base-base H-bond        (no_strong_hbond)
    stacking     stacking signature                                   (stacking)
    low_score    below the selection score floor                      (low_score)

`invalid` (failed geometric validation) and `adjacent` (sequence neighbours) are left
out, as they are in Figure 2A. As a consistency check the `selected` set of each
structure is compared with parse_pairs.json (the published output every figure uses).
Where they differ -- 6 of 1,183 structures on 2026-10-07: one pair each in 6QUL and
7K16 (engine build), and residue numbering/chains in the other four (a different
copy of the structure file) -- entries are listed under "fallback_to_output" and
retain selected output only in Figure 2/S3 agreement comparisons.
The script refuses to write if more than 1% differ.

Usage: python scripts/validation/extract_parse_candidates.py [PARSE_BINARY] [STRUCTURE_DIR]
Defaults: $PARSE_ENGINE or `parse` on PATH; data/structures (scripts/fetch_structures.py)
Writes data/comparison/parse_candidates.json:
    {pdb: {"status": "ok", "pairs": [[chain1, num1, chain2, num2, lw, disposition], ...]}}
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUTPUT_PAIRS = ROOT / "data" / "comparison" / "parse_pairs.json"
OUT = ROOT / "data" / "comparison" / "parse_candidates.json"
DEFAULT_BINARY = Path(os.environ.get("PARSE_ENGINE") or shutil.which("parse") or "parse")
DEFAULT_STRUCTURES = ROOT / "data" / "structures"
IDENTIFIED = {"selected", "competition", "ineligible", "stacking", "low_score"}


def split_res(res_id):
    """'R-4SU-8' -> ('R', '8'); the number keeps a sign and any insertion code."""
    chain, _name, num = res_id.split("-", 2)
    return chain, num


def candidates(binary, path):
    run = subprocess.run([str(binary), "dump-selection", str(path)],
                         capture_output=True, text=True, check=True)
    rows = []
    for line in run.stdout.splitlines():
        a, b, lw, disposition = line.split("\t")
        if disposition not in IDENTIFIED:
            continue
        (c1, n1), (c2, n2) = split_res(a), split_res(b)
        rows.append([c1, n1, c2, n2, None if lw == "None" else lw, disposition])
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("binary", nargs="?", type=Path, default=DEFAULT_BINARY)
    ap.add_argument("structures", nargs="?", type=Path, default=DEFAULT_STRUCTURES)
    ap.add_argument("--output", type=Path, default=OUT)
    args = ap.parse_args()
    binary, structures = args.binary.resolve(), args.structures.resolve()
    published = json.loads(OUTPUT_PAIRS.read_text())
    pdbs = sorted(p for p, r in published.items() if r.get("status") == "ok")

    def one(pdb):
        return pdb, candidates(binary, structures / f"{pdb}.cif")

    with ThreadPoolExecutor(max_workers=8) as pool:
        result = dict(pool.map(one, pdbs))

    def keyed_label(row):
        c1, n1, c2, n2, lw = row[:5]
        a, b = (c1, str(n1)), (c2, str(n2))
        labels = (lw or "?").split("|")
        if b < a:
            a, b = b, a
            labels = [x[0] + x[2] + x[1] if len(x) == 3 else x for x in labels]
        return (a, b), frozenset(labels)

    mismatched = []
    for pdb in pdbs:
        selected = dict(keyed_label(r) for r in result[pdb] if r[5] == "selected")
        output = dict(keyed_label(r) for r in published[pdb]["pairs"])
        if selected != output:
            changed = sum(selected.get(k) != output.get(k) for k in selected.keys() | output.keys())
            mismatched.append((pdb, changed))
    if len(mismatched) > 0.01 * len(pdbs):
        sys.exit(f"selected != parse_pairs.json for {len(mismatched)} structures, "
                 f"e.g. {mismatched[:5]}: the binary does not match the published output; "
                 "check both coordinates and dump-selection LW orientation")
    fallback = sorted(p for p, _ in mismatched)
    for pdb in fallback:   # published output only, all as "selected"
        result[pdb] = [[*r[:5], "selected"] for r in published[pdb]["pairs"]]

    data = {p: {"status": "ok", "pairs": result[p]} for p in pdbs}
    data["_meta"] = {
        "fallback_to_output": fallback,
        "identified_dispositions": sorted(IDENTIFIED),
        "label_order": "input_residue_order",
        "engine_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "note": "Use an engine whose dump-selection preserves residue/LW orientation. "
                "Figure 2/S3 use selected-output fallback; "
                "unclassified candidates are not identified-and-classified pairs.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, separators=(",", ":")))
    print(f"fallback to published output (no candidate credit): {fallback}")
    n = sum(len(v) for v in result.values())
    n_sel = sum(r[5] == "selected" for v in result.values() for r in v)
    print(f"wrote {args.output}: {len(pdbs)} structures, {n:,} "
          f"candidates ({n_sel:,} selected, {n - n_sel:,} dropped by selection)")


if __name__ == "__main__":
    main()
