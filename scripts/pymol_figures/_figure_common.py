"""Shared PyMOL helpers for the PARSE manuscript figure panels.

Not run directly — each fig*.py loads it. Keeps the rendering style identical
across figures: nucleotide colouring, explicit Watson-Crick H-bonds (not
PyMOL's promiscuous polar-contact search), white background, publication ray
settings.

Colours follow scripts/validation/view_redo_before_after.py, which produced the
published panels:  G = red, A = orange, C = forest, U = blue.
"""

import os

from pymol import cmd

def structures_dir():
    """Where the vendored figure coordinates live.

    Resolved relative to this file so the scripts work on any checkout. These
    are the minimal extracts committed by vendor_structures.py — the two
    residues each panel shows, plus 1GID chain A in full for panel 6A (~305 KB
    in total), so rebuilding a figure needs no download.
    """
    here = os.path.dirname(os.path.abspath(_FIGURE_COMMON_PATH))
    return os.path.normpath(
        os.path.join(here, "..", "..", "data", "figures", "structures"))

NT_RGB = {"G": [0.878, 0.0, 0.0], "A": [0.941, 0.502, 0.0],
          "C": [0.0, 0.565, 0.125], "U": [0.0, 0.0, 0.878]}

# Canonical Watson-Crick donor/acceptor atom pairs; first atom on residue 1.
WC_HBONDS = {
    ("G", "C"): [("N1", "N3"), ("N2", "O2"), ("O6", "N4")],
    ("C", "G"): [("N3", "N1"), ("O2", "N2"), ("N4", "O6")],
    ("A", "U"): [("N1", "N3"), ("N6", "O4")],
    ("U", "A"): [("N3", "N1"), ("O4", "N6")],
    ("G", "U"): [("O6", "N3"), ("N1", "O2")],
    ("U", "G"): [("N3", "O6"), ("O2", "N1")],
}
HB_CUTOFF = 3.6


def init_style():
    """Publication defaults: white background, thick sticks, ray tracing on."""
    for nt, rgb in NT_RGB.items():
        cmd.set_color("nt_" + nt, rgb)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("stick_radius", 0.18)
    cmd.set("dash_width", 2.5)
    cmd.set("dash_gap", 0.35)
    cmd.set("label_size", 16)
    cmd.set("label_color", "black")
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 0)
    cmd.set("depth_cue", 0)
    cmd.set("ray_shadows", 0)


def load_vendored(name, obj):
    """Load a vendored extract from data/figures/structures by file stem.

    Returns the object name, or None when the file is missing — the caller
    reports that rather than rendering an empty panel.
    """
    if obj in cmd.get_names("all"):
        return obj
    path = os.path.join(structures_dir(), name + ".pdb")
    if not os.path.exists(path):
        print("  !! missing vendored structure: %s" % path)
        print("     regenerate it with:  pymol -cq vendor_structures.py")
        return None
    cmd.load(path, obj)
    return obj


def actual_chain(obj, resi, want):
    """The chain id `resi` really has in `obj`.

    PDB format holds a single character, so a multi-character mmCIF chain (BA,
    iN) can be truncated when an extract is written. The vendored extracts hold
    only the pair's two residues, so resolving the chain from the file is both
    safe and immune to that.
    """
    found = set()
    cmd.iterate("%s and resi %s and name C1'" % (obj, resi),
                "found.add(chain)", space={"found": found})
    if want in found:
        return want
    if len(found) == 1:
        only = found.pop()
        print("      (chain %r not found; using %r from the file)" % (want, only))
        return only
    return want


def resn_of(obj, chain, resi):
    """Single-letter residue name of one residue, or '' if absent."""
    space = {"names": []}
    cmd.iterate("%s and chain %s and resi %s and name C1'" % (obj, chain, resi),
                "names.append(resn)", space=space)
    return space["names"][0].strip() if space["names"] else ""


def show_pair(obj, chain1, resi1, chain2, resi2, name, sticks_only=True):
    """Isolate one base pair as its own object, coloured by nucleotide.

    Returns the new object's name, or None when either residue is missing.
    """
    sel = ("(%s and chain %s and resi %s) or (%s and chain %s and resi %s)"
           % (obj, chain1, resi1, obj, chain2, resi2))
    if cmd.count_atoms(sel) == 0:
        print("  !! no atoms for %s %s/%s + %s/%s"
              % (obj, chain1, resi1, chain2, resi2))
        return None
    cmd.create(name, sel)
    cmd.hide("everything", name)
    cmd.show("sticks", name)
    if sticks_only:
        cmd.remove("%s and hydro" % name)
    for nt, _ in NT_RGB.items():
        cmd.color("nt_" + nt, "%s and resn %s+D%s" % (name, nt, nt))
    cmd.set("cartoon_side_chain_helper", 1, name)
    return name


def draw_wc_hbonds(pair_obj, chain1, resi1, chain2, resi2, label):
    """Draw only the canonical WC donor/acceptor contacts, when short enough.

    Deliberately not `cmd.distance(..., mode=2)`: PyMOL's polar-contact search
    invents contacts that are not the Watson-Crick set being illustrated.
    """
    n1 = resn_of(pair_obj, chain1, resi1)
    n2 = resn_of(pair_obj, chain2, resi2)
    pairs = WC_HBONDS.get((n1[-1:], n2[-1:]), [])
    drawn = 0
    for a1, a2 in pairs:
        s1 = "%s and chain %s and resi %s and name %s" % (pair_obj, chain1, resi1, a1)
        s2 = "%s and chain %s and resi %s and name %s" % (pair_obj, chain2, resi2, a2)
        if cmd.count_atoms(s1) == 0 or cmd.count_atoms(s2) == 0:
            continue
        d = cmd.get_distance(s1, s2)
        if d <= HB_CUTOFF:
            cmd.distance("%s_hb%d" % (label, drawn), s1, s2)
            drawn += 1
    for i in range(drawn):
        nm = "%s_hb%d" % (label, i)
        cmd.color("black", nm)
        cmd.hide("labels", nm)
    if drawn == 0:
        print("  (no canonical WC contact <= %.1f A — falling back to observed "
              "base-base polar contacts)" % HB_CUTOFF)
        return draw_base_polar_contacts(pair_obj, chain1, resi1,
                                        chain2, resi2, label)
    return drawn


# Base ring + exocyclic N/O atoms — the ones that can carry a base-base H-bond.
BASE_POLAR = ("N1+N2+N3+N4+N6+N7+N9+O2+O4+O6")


def draw_base_polar_contacts(pair_obj, chain1, resi1, chain2, resi2, label,
                             cutoff=3.5):
    """Draw observed base-base N/O contacts between the two residues.

    Needed for non-Watson-Crick classes (cWS, cWH, cSH ...), where the
    canonical WC atom table does not apply. Restricted to base polar atoms so
    the result is the base-base interaction being illustrated, not every
    backbone contact PyMOL can find.
    """
    s1 = ("%s and chain %s and resi %s and name %s"
          % (pair_obj, chain1, resi1, BASE_POLAR))
    s2 = ("%s and chain %s and resi %s and name %s"
          % (pair_obj, chain2, resi2, BASE_POLAR))
    if cmd.count_atoms(s1) == 0 or cmd.count_atoms(s2) == 0:
        return 0
    name = "%s_hb0" % label
    # cmd.distance returns the MEAN distance of the contacts it drew, not a
    # count, so report it as such.
    mean_d = cmd.distance(name, s1, s2, cutoff, mode=0)
    cmd.color("black", name)
    cmd.hide("labels", name)
    print("  observed base-base polar contacts <= %.1f A (mean %.2f A)"
          % (cutoff, mean_d))
    return 1


def frame(obj, buffer=3.0):
    cmd.orient(obj)
    cmd.zoom(obj, buffer)


def announce(panels):
    print("\n" + "=" * 62)
    print("panels: " + ", ".join(sorted(panels)))
    print('  panel("<name>")   set up one panel (hides the others)')
    print("=" * 62 + "\n")
