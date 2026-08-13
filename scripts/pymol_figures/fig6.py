"""Figure 6 panels A-C: PARSE quality scores mapped onto 1GID.

    (A) chain A coloured by base-pair quality (gray -> yellow -> orange -> red)
    (B) A235-U239 superposed on its idealized template   (PARSE score 51)
    (C) A C170 with B G254, likewise                     (PARSE score 91)

This is a LAUNCHER, not an implementation: every panel is produced by the
shipped PyMOL integration (parse_pymol.py, vendored alongside this file from the
engine's integrations/pymol/). Figure 6's claim is that PARSE launches PyMOL and
maps quality onto the structure, so the figure has to be made by that plugin and
not by a lookalike. All this file does is load the structure, load the plugin,
and name the panels so they match fig1 / figS1 / fig5.

Both chains are loaded because panel C's pair is BETWEEN the two copies of the
intron: 1GID crystallises as a dimer and C170 pairs with G254 of the other
chain, so a chain-A-only load cannot find it. (The manuscript caption calls it
"chain A of 1GID"; the contact is in fact A/C170-B/G254.)

Run me:
    python fig6.py            # opens PyMOL with panel 6A ready
    panel("6B")               # then File > Save Session As...

Needs the compiled engine: set PARSE_BINARY, or leave it and this script will
look for it in the sibling checkout.
"""

import os
import shutil
import sys

# Dual-mode entry point: `python fig6.py` launches PyMOL on this script (fresh
# session each run, so no `reinitialize` is needed); `run fig6.py` inside PyMOL
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
    """Directory holding this script (PyMOL's `run` does not set __file__)."""
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


HERE = _here()
PLUGIN = os.path.join(HERE, "parse_pymol.py")
# The plugin looks for idealized templates next to itself in the engine repo;
# this is a vendored copy, so point it at the vendored template set.
os.environ.setdefault("PARSE_IDEALS_DIR", os.path.normpath(os.path.join(
    HERE, "..", "..", "data", "figures", "basepair-idealized")))
STRUCT = os.path.normpath(os.path.join(
    HERE, "..", "..", "data", "figures", "structures", "1GID_chainAB.pdb"))
OBJ = "1GID_chainAB"

# panel -> (plugin command to run, what it shows)
PANELS = {
    "6A": (None, "chain A coloured by PARSE quality score"),
    "6B": ("A-235", "A235-U239 vs idealized template (score 51)"),
    "6C": ("A-170", "A/C170 with B/G254 vs idealized template (score 91)"),
}


def _find_binary():
    """Locate the compiled engine, preferring an explicit PARSE_BINARY."""
    if os.environ.get("PARSE_BINARY"):
        return os.environ["PARSE_BINARY"]
    guess = os.path.normpath(os.path.join(
        HERE, "..", "..", "..", "pair-finder-main-cpp", "build", "parse"))
    if os.path.exists(guess):
        return guess
    return shutil.which("parse")


def build():
    """Load the structure, load the plugin, and score the whole thing."""
    cmd.delete("all")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)

    if not os.path.exists(STRUCT):
        print("!! missing %s — run: pymol -cq vendor_structures.py" % STRUCT)
        return
    cmd.load(STRUCT, OBJ)

    binary = _find_binary()
    if binary:
        os.environ["PARSE_BINARY"] = binary
        print("engine: %s" % binary)
    else:
        print("!! `parse` binary not found — set PARSE_BINARY to the compiled "
              "engine, then run build() again")
        return

    if "parse_score" not in cmd.keyword:
        cmd.do("run " + PLUGIN)          # no quotes: run takes a bare path
    if "parse_score" not in cmd.keyword:
        print("!! could not load the plugin at %s" % PLUGIN)
        return

    # tier=all so both Review and Acceptable pairs carry colour, matching the
    # caption's gray -> yellow -> orange -> red gradient.
    cmd.do("parse_score %s, all" % OBJ)
    print("\n" + "=" * 62)
    print('  panel("6A") / panel("6B") / panel("6C")')
    print("  then save the session yourself:  File > Save Session As...")
    print("=" * 62 + "\n")


def panel(name="6A"):
    """Set up one Figure-6 panel by delegating to the plugin."""
    if OBJ not in cmd.get_names("all"):
        build()
    if name not in PANELS:
        print("unknown panel %r; have 6A, 6B, 6C" % name)
        return
    goto, desc = PANELS[name]
    if "parse_ideal" not in cmd.keyword:
        print("!! plugin not loaded — run build() first")
        return
    if goto is None:
        cmd.do("parse_ideal off")
        cmd.enable(OBJ)
        # Hide chain B EVERYWHERE, not just in the structure: parse_score also
        # builds `parse_overlay`, which carries the flagged pairs of both
        # chains, so hiding it only in OBJ leaves the second copy's sticks
        # floating in the view.
        cmd.hide("everything", "chain B")
        cmd.show("cartoon", "%s and chain A" % OBJ)
        cmd.orient("%s and chain A" % OBJ)
        cmd.zoom("%s and chain A" % OBJ, 2)
    else:
        cmd.show("cartoon", "%s and chain B" % OBJ)   # panel C spans both copies
        cmd.do("parse_ideal on")
        cmd.do("parse_goto %s" % goto)
        # The figure shows the pair alone with its ideal superposed, so drop the
        # whole structure, the all-pairs overlay and the surrounding context —
        # leaving only parse_focus (the pair) and parse_ideal_pair (the ghost).
        cmd.do("disable %s" % OBJ)
        cmd.do("disable parse_overlay")
        cmd.do("disable parse_context")
        # The plugin annotates the pair for interactive use: an atom label on
        # parse_focus and a viewport title. Neither belongs in the figure.
        cmd.hide("labels", "parse_focus")
        cmd.hide("labels", "parse_ideal_pair")
        try:
            cmd.set_title(OBJ, 1, "")
        except Exception:  # noqa: BLE001 - cosmetic only
            pass
        cmd.orient("parse_focus")
        cmd.zoom("parse_focus", 2)
    print("%s: %s" % (name, desc))


def ideal(state="toggle"):
    """Show/hide the green idealized-template ghost on the current panel.

    Lets you save the same pair with and without its reference geometry:
        ideal("off")   hide it      ideal("on")   show it      ideal()  toggle
    """
    objs = [o for o in ("parse_ideal_pair", "parse_ideal_hb")
            if o in cmd.get_names("objects")]
    if not objs:
        print("no idealized template loaded — run panel(\"6B\") or panel(\"6C\")")
        return
    want = str(state).lower()
    if want in ("toggle", ""):
        want = "off" if cmd.get_object_settings("parse_ideal_pair") is not None \
            and "parse_ideal_pair" in cmd.get_names("objects", enabled_only=1) else "on"
    for o in objs:
        (cmd.enable if want == "on" else cmd.disable)(o)
    print("idealized template %s" % ("shown" if want == "on" else "hidden"))


cmd.extend("ideal", ideal)
cmd.extend("panel", panel)
cmd.extend("build", build)
build()
panel("6A")
