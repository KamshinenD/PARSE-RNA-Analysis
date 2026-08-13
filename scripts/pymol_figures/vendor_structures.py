"""Vendor the minimal coordinates the PyMOL figure panels need.

Run ONCE (by us, not by readers) to populate
``data/figures/structures/``. After that the figure scripts load from the repo
and need no network.

Why extracts rather than whole entries: 7A0S is 13 MB and 6CAQ 8.7 MB, but every
panel except 6A shows exactly two residues. Extracting them keeps the vendored
set at a few hundred KB, which is reasonable to commit; whole ribosomes are not.
1GID is kept as a full chain because panel 6A scores the entire chain.

    pymol -cq vendor_structures.py
"""

import os

from pymol import cmd

REPO = ("/Users/kdewan2/Desktop/Projects/find-pair-score-rna/PARSE-RNA-Analysis")
OUT = os.path.join(REPO, "data", "figures", "structures")
LOCAL_PDB = "/Users/kdewan2/Desktop/Projects/find_pair_2/data/pdb"
REDO_DIR = ("/Users/kdewan2/Desktop/Projects/find-pair-score-rna/"
            "prototyped-pair-finder-main/data/validation/redo_files")

# name -> (pdb id, selection, use the PDB-REDO copy?)
TARGETS = [
    # Figure 6: whole chain A, because panel 6A scores it.
    ("1GID_chainA",        "1GID", "chain A", False),
    # Figure 5C/D: one pair each, original and re-refined.
    ("4V2S_Q11_Q41_orig",  "4V2S", "chain Q and resi 11+41", False),
    ("4V2S_Q11_Q41_redo",  "4V2S", "chain Q and resi 11+41", True),
    ("6ZDP_B1296_B1322_orig", "6ZDP", "chain B and resi 1296+1322", False),
    ("6ZDP_B1296_B1322_redo", "6ZDP", "chain B and resi 1296+1322", True),
    # Supplementary S1B/C/D: one pair each.
    ("7A0S_X1100_X1113",   "7A0S", "chain X and resi 1100+1113", False),
    ("6CAQ_A1036_A1004",   "6CAQ", "chain A and resi 1036+1004", False),
    ("5J8B_A2713_A2715",   "5J8B", "chain A and resi 2713+2715", False),
]


def source_path(pdb_id, redo):
    p = (os.path.join(REDO_DIR, pdb_id, pdb_id + ".pdb") if redo
         else os.path.join(LOCAL_PDB, pdb_id + ".pdb"))
    return p if os.path.exists(p) else None


def main():
    os.makedirs(OUT, exist_ok=True)
    total = 0
    for name, pdb_id, sele, redo in TARGETS:
        cmd.delete("all")
        src = source_path(pdb_id, redo)
        if src:
            cmd.load(src, "src")
        else:
            if redo:
                print("!! %s: no local PDB-REDO copy; skipping" % name)
                continue
            print("   %s: not local, fetching from RCSB" % pdb_id)
            cmd.fetch(pdb_id, "src", async_=0)
        if cmd.count_atoms("src") == 0:
            print("!! %s: nothing loaded" % name)
            continue
        n = cmd.count_atoms("src and (%s)" % sele)
        if n == 0:
            print("!! %s: selection '%s' matched no atoms" % (name, sele))
            continue
        cmd.create("out", "src and (%s)" % sele)
        cmd.remove("out and solvent")
        path = os.path.join(OUT, name + ".pdb")
        cmd.save(path, "out")
        kb = os.path.getsize(path) / 1024.0
        total += kb
        print("   %-26s %6d atoms  %7.1f KB" % (name, n, kb))
    print("\nvendored %.1f KB total -> %s" % (total, OUT))


main()
