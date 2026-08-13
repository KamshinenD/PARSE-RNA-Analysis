"""Supplementary Figure S1 panels B-D: the three kinds of PARSE/DSSR disagreement.

Panel A is the stacked bar chart of disagreement composition. B, C and D give one
worked structural example of each of its three categories, so a reader can see
what each bar actually looks like:

    (B) Different LW assignment     — both tools find the pair but call a
                                      different interaction edge.
                                      4DR6 chain A, A1055-U1205
                                      PARSE cWS / DSSR cWW, coplanar to 0.5 deg
                                      with a single H-bond
    (C) Detected only by PARSE      — a real pair DSSR does not report.
                                      3CPW chain 0, C2071-U625: cWW,
                                      PARSE score 86.7
    (D) Excluded during candidate   — DSSR calls it a pair; the bases are 9 A
        generation                    apart, 59 deg twisted and share no
                                      H-bond, so PARSE never scores it.
                                      6YAL chain 2, C1585-G1647

Examples are drawn from data/comparison/<class>.json, whose `reason` field is
exactly what panel A is binned on: `dssr_assigns:<class>` -> B, `novel` -> C,
`not_candidate` -> D.

Run me:
    python figS1.py           # opens PyMOL with panel S1B ready
    panel("S1C")
    panel("all")              # all three at once
"""

import os
import shutil
import sys

# Dual-mode entry point: `python figS1.py` launches PyMOL on this script (fresh
# session each run, so no `reinitialize` is needed); `run figS1.py` inside PyMOL
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

PANELS = {
    "S1B": {
        "file": "figS1B_edge_4DR6", "pdb": "4DR6",
        "chain1": "A", "resi1": "1055", "chain2": "A", "resi2": "1205",
        "category": "Different LW assignment",
        "detail": "PARSE cWS vs DSSR cWW — coplanar (0.5 deg), one H-bond, "
                  "so which edge is engaged is genuinely ambiguous",
        "hbonds": True,
    },
    "S1C": {
        "file": "figS1C_novel_3CPW", "pdb": "3CPW",
        "chain1": "0", "resi1": "2071", "chain2": "0", "resi2": "625",
        "category": "Detected only by PARSE",
        "detail": "cWW C-U, PARSE score 86.7 — a pair DSSR does not report",
        "hbonds": True,
    },
    "S1D": {
        "file": "figS1D_rejected_6YAL", "pdb": "6YAL",
        "chain1": "2", "resi1": "1585", "chain2": "2", "resi2": "1647",
        "category": "Excluded during candidate generation",
        "detail": "DSSR calls this cSS; the bases are 9 A apart and 59 deg "
                  "twisted, with no H-bond, so PARSE never scores it",
        # Drawing a bond here would assert the very thing the panel shows is
        # absent; the geometry is printed instead.
        "hbonds": False,
    },
}

_built = []


def _report_geometry(obj, ch1, r1, ch2, r2):
    """Print the base-plane geometry that made PARSE reject this pair.

    Panel D is a pair DSSR reports and PARSE does not, so the useful caption
    number is the geometry, not a bond distance: the two bases sit far apart
    and twisted, which is what the candidate stage screens on.
    """
    import numpy as np

    purine = ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"]
    pyrim = ["N1", "C2", "N3", "C4", "C5", "C6"]

    def plane(ch, resi):
        names = set()
        cmd.iterate("%s and chain %s and resi %s" % (obj, ch, resi),
                    "names.add(name)", space={"names": names})
        ring = purine if "N9" in names else pyrim
        pts = []
        for n in ring:
            m = cmd.get_model("%s and chain %s and resi %s and name %s"
                              % (obj, ch, resi, n))
            if m.atom:
                pts.append(m.atom[0].coord)
        if len(pts) < 3:
            return None, None
        P = np.array(pts)
        c = P.mean(axis=0)
        return c, np.linalg.svd(P - c)[2][2]

    c1, n1 = plane(ch1, r1)
    c2, n2 = plane(ch2, r2)
    if c1 is None or c2 is None:
        return
    ang = np.degrees(np.arccos(abs(float(np.dot(n1, n2)))))
    sep = float(np.linalg.norm(c2 - c1))
    print("      base planes %.0f deg apart, ring centroids %.1f A apart "
          "— too far and too twisted to be a candidate pair" % (ang, sep))


def build():
    init_style()
    cmd.delete("all")
    _built[:] = []
    for name in sorted(PANELS):
        c = PANELS[name]
        print("[%s] %s" % (name, c["category"]))
        print("      %s  |  %s %s/%s + %s/%s"
              % (c["detail"], c["pdb"], c["chain1"], c["resi1"],
                 c["chain2"], c["resi2"]))
        src = load_vendored(c["file"], "src_" + name)
        if src is None:
            continue
        ch1 = actual_chain(src, c["resi1"], c["chain1"])
        ch2 = actual_chain(src, c["resi2"], c["chain2"])
        obj = show_pair(src, ch1, c["resi1"], ch2, c["resi2"], name)
        if obj is None:
            cmd.delete(src)
            continue
        if c["hbonds"]:
            draw_wc_hbonds(obj, ch1, c["resi1"], ch2, c["resi2"], name)
        else:
            _report_geometry(obj, ch1, c["resi1"], ch2, c["resi2"])
        cmd.delete(src)
        _built.append(name)
    announce(_built)


def panel(name="S1B"):
    """Show one panel (or 'all'), hiding the others, and frame it."""
    if not _built:
        build()
    for n in _built:
        cmd.disable(n)
        for i in range(4):
            cmd.disable("%s_hb%d" % (n, i))
    if name == "all":
        for n in _built:
            cmd.enable(n)
            for i in range(4):
                cmd.enable("%s_hb%d" % (n, i))
        cmd.zoom("all", 3)
        return
    if name not in _built:
        print("unknown panel %r; have: %s" % (name, ", ".join(_built)))
        return
    cmd.enable(name)
    for i in range(4):
        cmd.enable("%s_hb%d" % (name, i))
    frame(name)
    c = PANELS[name]
    print("%s: %s — %s" % (name, c["category"], c["detail"]))


cmd.extend("panel", panel)
cmd.extend("build", build)
build()
panel("S1B")
