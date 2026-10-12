"""Base identity of every residue in the tool-comparison structures.

Figure 2 splits base-pair detection agreement into canonical (Watson-Crick G-C / A-U)
and non-canonical pairs. The tool pair lists in data/comparison/*_pairs.json carry
only residue ids and LW labels, and the residue names that are on disk (FastRNA's
pairs.csv, FR3D's unit ids) cover only the residues those tools pair -- 8-26% of
the other tools' pairs would go unclassified. This reads the structures themselves
once and writes the parent base (A/C/G/U, modified residues mapped to their parent
by gemmi's residue table) of every model-1 residue that has one, keyed the way the
figure notebooks key residues: "chain|number|icode". Residues without an A/C/G/U
parent (protein, ligands, unknown) are left out, so a missing key means "not a
standard base" and such a pair counts as non-canonical.

Usage: python scripts/validation/extract_residue_names.py [STRUCTURE_DIR]
Default STRUCTURE_DIR: data/structures (fetch with scripts/fetch_structures.py)
Writes data/comparison/residue_names.json.
"""
import json
import sys
from pathlib import Path

import gemmi

ROOT = Path(__file__).resolve().parents[2]
PAIRS = ROOT / "data" / "comparison" / "parse_pairs.json"
OUT = ROOT / "data" / "comparison" / "residue_names.json"
DEFAULT_STRUCTURES = ROOT / "data" / "structures"


def parent_base(name):
    if name in ("A", "C", "G", "U"):
        return name
    info = gemmi.find_tabulated_residue(name)
    code = info.one_letter_code.upper() if info is not None else ""
    return code if code in ("A", "C", "G", "U") else "?"


def residue_names(path):
    st = gemmi.read_structure(str(path))
    out = {}
    for chain in st[0]:
        for res in chain:
            base = parent_base(res.name)
            if base != "?":
                out[f"{chain.name}|{res.seqid.num}|{res.seqid.icode.strip()}"] = base
    return out


def main():
    structures = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_STRUCTURES
    pdb_ids = sorted(json.loads(PAIRS.read_text()))
    result, missing = {}, []
    for i, pdb in enumerate(pdb_ids, 1):
        path = next((p for p in (structures / f"{pdb}.cif", structures / f"{pdb.lower()}.cif")
                     if p.exists()), None)
        if path is None:
            missing.append(pdb)
            continue
        result[pdb] = residue_names(path)
        if i % 100 == 0:
            print(f"{i}/{len(pdb_ids)}", flush=True)
    if missing:
        sys.exit(f"no structure file for {len(missing)} entries, e.g. {missing[:5]}")
    OUT.write_text(json.dumps(result, separators=(",", ":")))
    print(f"wrote {OUT.relative_to(ROOT)}: {len(result)} structures, "
          f"{sum(map(len, result.values())):,} residues")


if __name__ == "__main__":
    main()
