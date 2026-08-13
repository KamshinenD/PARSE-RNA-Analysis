# PyMOL figure panels

Each script opens its own PyMOL session with the figure's panels prepared. Pick a
panel, orient it, then save the session yourself (**File > Save Session As…**).

```bash
cd scripts/pymol_figures
python fig1.py        # Figure 1  A, C
python figS1.py       # Supplementary S1  B, C, D
python fig5.py        # Figure 5  C, D
```

Figure 6 is produced by the engine's PyMOL plugin rather than a script of its
own — see below.

All coordinates are vendored in `../../data/figures/structures/` (minimal
extracts, ~1 MB total), so nothing is downloaded at run time.

---

## Figure 1 — ideal vs distorted base-pair geometry

```
python fig1.py
panel("1A")                 # canonical cWW G-C
panel("1C_ideal_GC")        # 9AXU 2/G2456-C2467   ProSco 100 on all six
panel("1C_ideal_AU")        # 7UNW A/A42-U428      five at 100, opening 92
panel("1C_distorted_GC")    # 1NJP 0/G726-C730     all ~0, 1 H-bond
panel("1C_distorted_AU")    # 7D6Z f/A81-U88       propeller -65.7 deg, 1 H-bond
```

"Ideal" is a real pair selected by **minimum ProSco across all six geometry
parameters**, not by score 100 — score 100 only needs ProSco >= 5 and still
admits visible distortion. Real pairs are used rather than the idealized
templates because those hold only base atoms + C1' (21 atoms, no sugar or
phosphate).

## Supplementary Figure S1 — the three kinds of PARSE/DSSR disagreement

```
python figS1.py
panel("S1B")   # 4DR6 A/A1055-U1205  PARSE cWS vs DSSR cWW, coplanar 0.5 deg, 1 H-bond
panel("S1C")   # 3CPW 0/C2071-U625   cWW, score 86.7, not reported by DSSR
panel("S1D")   # 6YAL 2/C1585-G1647  DSSR cSS; 9 A apart, 59 deg twisted, no H-bond
panel("all")   # all three at once
```

Panels map to the three categories of panel A's bar chart, taken from
`data/comparison/<class>.json` (`dssr_assigns:` -> B, `novel` -> C,
`not_candidate` -> D).

## Figure 5 — before/after PDB-REDO

```
python fig5.py
panel("5C")           # 4V2S Q/G11-C41      original
panel("5C", "redo")   # the same pair after re-refinement
panel("5D")           # 6ZDP B/G1296-C1322  original
panel("5D", "redo")
```

The two copies are separate objects shown one at a time, never superimposed.

## Figure 6 — score visualisation, via the engine plugin

Figure 6 *is* the shipped PyMOL integration, so it is generated with the plugin
itself (`parse_pymol.py`, vendored here from the engine's
`integrations/pymol/`) rather than a wrapper script:

```bash
export PARSE_BINARY=/path/to/PARSE-RNA/build/parse
pymol ../../data/figures/structures/1GID_chainAB.pdb
```

then inside PyMOL:

```
run parse_pymol.py
parse_score 1GID_chainAB, all      # panel A — colour by quality score
parse_ideal on
parse_goto A-235                   # panel B — A235-U239, score 51
parse_goto A-170                   # panel C — A C170 with B G254, score 91
```

Panel A shows chain A in the manuscript; chain B is loaded because **panel C's
pair is between the two copies of the intron** — 1GID crystallises as a dimer
and C170 pairs with G254 of the *other* chain, so a chain-A-only load cannot
find it. (The manuscript caption calls it "chain A of 1GID"; the contact is
in fact A/C170-B/G254.)

Useful plugin commands: `parse_list`, `parse_next` / `parse_prev`,
`parse_info`, `parse_clear`, `parse_ideal off`.

---

## Regenerating the vendored coordinates

```bash
pymol -cq vendor_structures.py
```

Extracts only what the panels need — two residues per pair, plus 1GID chains A
and B for Figure 6.
