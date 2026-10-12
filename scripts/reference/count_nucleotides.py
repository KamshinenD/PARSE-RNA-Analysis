"""Count nucleotides per structure and write them into the PDB metadata cache.

    python scripts/reference/count_nucleotides.py            # update metadata_cache.json
    python scripts/reference/count_nucleotides.py --dry-run  # report changes only

A nucleotide is any residue of model 1 that belongs to a polymer entity the entry
declares as RNA, DNA or an RNA/DNA hybrid (mmCIF `_entity_poly.type`). The residue
name is not consulted, so modified (e.g. pseudouridine, 2'-O-methyl), locked (LNA)
and mirror-image (L-RNA) nucleotides are all counted, and bound nucleotide ligands
(non-polymer entities) are not.

This replaces the original count, which kept only residues named A/C/G/U, DA/DC/DG/DT
and eight common modified bases: structures built entirely from other residue types
(e.g. 2X2Q, LNA; 8F24, L-RNA; 469D, 2'-O-methoxyethyl RNA) were counted as 0 and
modified residues elsewhere were left out.

Reads data/structures/<PDB>.cif (scripts/fetch_structures.py); entries without a local
file keep their existing count. Updates nucleotide_counts in
data/reference/metadata/metadata_cache.json in place.
"""
from __future__ import annotations

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import gemmi

REPO = Path(__file__).resolve().parents[2]
STRUCTURES = REPO / "data" / "structures"
METADATA = REPO / "data" / "reference" / "metadata" / "metadata_cache.json"
NUCLEIC = {gemmi.PolymerType.Rna, gemmi.PolymerType.Dna, gemmi.PolymerType.DnaRnaHybrid}


def count(pdb: str) -> tuple[str, int | None]:
    path = STRUCTURES / f"{pdb}.cif"
    if not path.exists():
        return pdb, None
    try:
        st = gemmi.read_structure(str(path))
        st.setup_entities()
        seen = set()
        for chain in st[0]:
            for res in chain:
                if res.entity_type != gemmi.EntityType.Polymer or not res.subchain:
                    continue
                entity = st.get_entity_of(chain.get_subchain(res.subchain))
                if entity is not None and entity.polymer_type in NUCLEIC:
                    seen.add((chain.name, str(res.seqid)))
        return pdb, len(seen)
    except Exception:          # unreadable file: keep the existing count
        return pdb, None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="report changes, write nothing")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    meta = json.loads(METADATA.read_text())
    counts = meta["nucleotide_counts"]
    with ProcessPoolExecutor(args.workers) as ex:
        results = dict(ex.map(count, list(counts), chunksize=20))

    changed = {p: (counts[p], n) for p, n in results.items() if n is not None and counts.get(p) != n}
    counted = sum(n is not None for n in results.values())
    print(f"{counted:,} structures counted from local files "
          f"({len(counts) - counted:,} without a file keep their count)")
    print(f"{len(changed):,} counts changed; zeros before {sum(counts[p] == 0 for p in results if results[p] is not None)}"
          f" -> after {sum(n == 0 for n in results.values() if n is not None)}")
    for p, (old, new) in sorted(changed.items(), key=lambda kv: kv[1][1] - kv[1][0])[-5:]:
        print(f"  {p}: {old} -> {new}")
    if args.dry_run:
        return
    for p, (_, new) in changed.items():
        counts[p] = new
    METADATA.write_text(json.dumps(meta))
    print(f"wrote {METADATA.relative_to(REPO)}")


if __name__ == "__main__":
    main()
