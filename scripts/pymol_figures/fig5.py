"""Figure 5 panels C-D: the same base pair before and after PDB-REDO refinement.

Identities from the claims-document caption:
    (C) 4V2S chain Q  G11-C41      PARSE 26.4 -> 100
    (D) 6ZDP chain B  G1296-C1322  PARSE 28.3 -> 100

The original and re-refined copies are kept as SEPARATE objects and shown one
at a time (never superimposed), matching
scripts/validation/view_redo_before_after.py, which produced the published
panels.

Usage (inside PyMOL):
    run /Users/kdewan2/Desktop/Projects/find-pair-score-rna/PARSE-RNA-Analysis/scripts/pymol_figures/fig5.py
    panel("5C")          # original copy
    panel("5C", "redo")  # re-refined copy
    panel("5D") ; panel("5D", "redo")

Panels A and B are plots (figure5.ipynb), not PyMOL.
"""

import os
import shutil
import sys

# Dual-mode entry point.
#   python figX.py     -> launches PyMOL on this script (fresh session each run,
#                         so no `reinitialize` is ever needed)
#   run figX.py        -> already inside PyMOL, carry on
# `pymol` is importable in this environment even from plain python, so the test
# is whether it is already LOADED, not whether it can be imported.
if "pymol" not in sys.modules:
    _self = os.path.abspath(__file__)
    # PyMOL's `run` does not set __file__, so pass our directory through the
    # environment for the relaunched process to pick up.
    os.environ["PARSE_FIG_DIR"] = os.path.dirname(_self)
    _exe = shutil.which("pymol")
    if _exe is None:
        sys.exit("PyMOL not found on PATH — activate the environment that has it.")
    os.execv(_exe, [_exe, _self])

from pymol import cmd


def _here():
    """Directory holding these figure scripts."""
    d = os.environ.get("PARSE_FIG_DIR")
    if d and os.path.isdir(d):
        return d
    try:
        import pymol
        p = getattr(pymol, "__script__", None)
        if p:
            return os.path.dirname(os.path.abspath(p))
    except Exception:  # noqa: BLE001
        pass
    return os.getcwd()


_FIGURE_COMMON_PATH = os.path.join(_here(), "_figure_common.py")
exec(open(_FIGURE_COMMON_PATH).read())


PANELS = {
    "5C": {"file": "4V2S_Q11_Q41", "pdb": "4V2S", "chain1": "Q", "resi1": "11",
           "chain2": "Q", "resi2": "41", "orig": 26.4},
    "5D": {"file": "6ZDP_B1296_B1322", "pdb": "6ZDP", "chain1": "B", "resi1": "1296",
           "chain2": "B", "resi2": "1322", "orig": 28.3},
}

_built = []


def build():
    """Build an `<panel>_orig` and `<panel>_redo` object for each case."""
    init_style()
    cmd.delete("all")
    _built[:] = []
    for name in sorted(PANELS):
        c = PANELS[name]
        print("[%s] %s %s/%s + %s/%s   score %.1f -> 100"
              % (name, c["pdb"], c["chain1"], c["resi1"],
                 c["chain2"], c["resi2"], c["orig"]))
        for tag, redo in (("orig", False), ("redo", True)):
            src = load_vendored("%s_%s" % (c["file"], tag), "src_%s_%s" % (name, tag))
            if src is None:
                print("  !! %s copy unavailable" % tag)
                continue
            obj = show_pair(src, c["chain1"], c["resi1"],
                            c["chain2"], c["resi2"], "%s_%s" % (name, tag))
            if obj is None:
                cmd.delete(src)
                continue
            n = draw_wc_hbonds(obj, c["chain1"], c["resi1"],
                               c["chain2"], c["resi2"], "%s_%s" % (name, tag))
            print("  %-4s : %d canonical WC H-bond(s) <= %.1f A" % (tag, n, HB_CUTOFF))
            cmd.delete(src)
            _built.append("%s_%s" % (name, tag))
    announce(_built)


def panel(name="5C", which="orig"):
    """Show `<name>_<which>` alone and frame it. which = orig | redo."""
    if not _built:
        build()
    target = "%s_%s" % (name, which)
    for n in _built:
        cmd.disable(n)
        for i in range(4):
            cmd.disable("%s_hb%d" % (n, i))
    if target not in _built:
        print("unknown panel %r; have: %s" % (target, ", ".join(_built)))
        return
    cmd.enable(target)
    for i in range(4):
        cmd.enable("%s_hb%d" % (target, i))
    frame(target)
    c = PANELS[name]
    label = "original (score %.1f)" % c["orig"] if which == "orig" else "PDB-REDO (score 100)"
    print("%s  %s  %s" % (name, c["pdb"], label))


cmd.extend("panel", panel)
cmd.extend("build", build)
build()
panel("5C", "orig")
