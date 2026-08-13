"""Figure 1 panels A and C: ideal vs severely distorted base-pair geometry.

    (A) a canonical cWW G-C pair
    (C) upper  — ideal cWW G-C and A-U
        lower  — severely distorted cWW G-C and A-U

Panel B (backbone torsions) is a schematic, not PyMOL, so it is not here.

Panels are selected by ProSco (the per-parameter density percentile), not by
the pair score: score 100 only requires ProSco >= 5, which still admits visible
distortion. The ideal G-C here is the best in the 368,169-pair reference set,
with ProSco 100 on all six geometry parameters. The distorted examples sit at
the bottom of every distribution AND have lost all but one base-base hydrogen
bond, so the damage is visible as missing H-bonds rather than a slight tilt.

Real pairs are used rather than the idealized templates because those hold only
base atoms + C1' (21 atoms, no ribose or phosphate); these panels need the
sugar and backbone.

Run me:
    python fig1.py            # opens PyMOL with panel 1A ready
    panel("1C_ideal_GC")
"""

import os
import shutil
import sys

# Dual-mode entry point: `python fig1.py` launches PyMOL on this script (fresh
# session each run, so no `reinitialize` is needed); `run fig1.py` inside PyMOL
# just carries on. pymol is importable from plain python here, so the test is
# whether it is already LOADED.
if "pymol" not in sys.modules:
    _self = os.path.abspath(__file__)
    os.environ["PARSE_FIG_DIR"] = os.path.dirname(_self)
    _exe = shutil.which("pymol")
    if _exe is None:
        sys.exit("PyMOL not found on PATH — activate the environment that has it.")
    os.execv(_exe, [_exe, _self])

from pymol import cmd


def _here():
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

# Panels are chosen by ProSco, the density percentile of each parameter — not
# merely by "score 100", which only requires ProSco >= 5 (the Preferred cut).
#
# PERFECT: ranked by the MINIMUM ProSco across all six geometry parameters, so
# the pair's *worst* parameter is still at the top of its distribution. The G-C
# is the single best in the 368k-pair reference set: 100 on every parameter.
PERFECT = {
    "1A":          {"file": "fig1_perfect_GC_9AXU", "chain1": "2",
                    "resi1": "2456", "chain2": "2", "resi2": "2467",
                    "seq": "G-C", "pdb": "9AXU",
                    "prosco": "all six = 100"},
    "1C_ideal_GC": {"file": "fig1_perfect_GC_9AXU", "chain1": "2",
                    "resi1": "2456", "chain2": "2", "resi2": "2467",
                    "seq": "G-C", "pdb": "9AXU",
                    "prosco": "all six = 100"},
    "1C_ideal_AU": {"file": "fig1_perfect_AU_7UNW", "chain1": "A",
                    "resi1": "42", "chain2": "A", "resi2": "428",
                    "seq": "A-U", "pdb": "7UNW",
                    "prosco": "five = 100, opening = 92 (best A-U available)"},
}

# DISTORTED: the opposite extreme — every parameter at the bottom of its
# distribution AND the H-bond complement collapsed to a single base-base bond,
# so the damage is visible as missing hydrogen bonds, not just a tilt.
DISTORTED = {
    "1C_distorted_GC": {"file": "fig1_distorted_GC_1NJP", "chain1": "0",
                        "resi1": "726", "chain2": "0", "resi2": "730",
                        "seq": "G-C", "pdb": "1NJP", "hbonds": 1,
                        "prosco": "shear 6, stretch 0, stagger 1, buckle 0, "
                                  "propeller 0, opening 0"},
    # Twisted out of plane rather than pulled apart: propeller -65.7 deg
    # (ProSco 0) with buckle -16.1 deg, so the two bases are rotated against
    # each other and only one H-bond survives — the geometry the published
    # panel shows.
    "1C_distorted_AU": {"file": "fig1_distorted_AU_7D6Z", "chain1": "f",
                        "resi1": "81", "chain2": "f", "resi2": "88",
                        "seq": "A-U", "pdb": "7D6Z", "hbonds": 1,
                        "prosco": "propeller -65.7 deg (ProSco 0), "
                                  "buckle -16.1 deg (6) — twisted"},
}

_built = []


def _build_pair(name, c):
    """One panel: a cWW pair from its own vendored extract."""
    src = load_vendored(c["file"], "src_" + name)
    if src is None:
        return
    ch1 = actual_chain(src, c["resi1"], c["chain1"])
    ch2 = actual_chain(src, c["resi2"], c["chain2"])
    obj = show_pair(src, ch1, c["resi1"], ch2, c["resi2"], name)
    if obj is None:
        return
    n = draw_wc_hbonds(obj, ch1, c["resi1"], ch2, c["resi2"], name)
    print("      %d canonical WC H-bond(s) drawn" % n)
    cmd.delete(src)
    _built.append(name)


def build():
    init_style()
    cmd.delete("all")
    _built[:] = []
    for name, c in PERFECT.items():
        print("[%s] %s %s/%s-%s  %s  ProSco: %s"
              % (name, c["pdb"], c["chain1"], c["resi1"], c["resi2"],
                 c["seq"], c["prosco"]))
        _build_pair(name, c)
    for name, c in DISTORTED.items():
        print("[%s] %s %s/%s-%s  %s  %d H-bond  ProSco: %s"
              % (name, c["pdb"], c["chain1"], c["resi1"], c["resi2"],
                 c["seq"], c["hbonds"], c["prosco"]))
        _build_pair(name, c)
    announce(_built)


def panel(name="1A"):
    """Show one panel alone and frame it."""
    if not _built:
        build()
    for n in _built:
        cmd.disable(n)
        for i in range(4):
            cmd.disable("%s_hb%d" % (n, i))
    if name not in _built:
        print("unknown panel %r; have: %s" % (name, ", ".join(_built)))
        return
    cmd.enable(name)
    for i in range(4):
        cmd.enable("%s_hb%d" % (name, i))
    frame(name)
    c = DISTORTED.get(name) or PERFECT[name]
    kind = "severely distorted" if name in DISTORTED else "ideal"
    print("%s: %s %s/%s-%s %s  (%s)  ProSco: %s"
          % (name, c["pdb"], c["chain1"], c["resi1"], c["resi2"], c["seq"],
             kind, c["prosco"]))


cmd.extend("panel", panel)
cmd.extend("build", build)
build()
panel("1A")
