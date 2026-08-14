"""Save every structural panel as its own .pse, into one folder.

    pymol -cq save_all_pses.py

Writes figures/pse/<panel>.pse for all panels of Figures 1, 5, 6 and
Supplementary S1, then you can zip that folder and send it.

Each file holds ONE panel. A PyMOL session normally stores every object it has,
including hidden ones, so a naive save of `1A` would also embed the other four
Figure-1 panels. Here the panel is set up first and everything still disabled is
deleted, so what lands in the file is only what is on screen.

The camera is whatever the script's default framing gives. Re-orienting by hand
and re-saving any individual panel is still worth doing for a final figure; this
is for handing over the whole set at once.
"""

import glob
import os

from pymol import cmd

HERE = os.path.dirname(os.path.abspath(
    __import__("pymol").__script__ if hasattr(__import__("pymol"), "__script__")
    else __file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "..", "figures", "pse"))

# script -> panels to save. fig5's panels take a second argument.
JOBS = [
    ("fig1.py", ["1A", "1C_ideal_GC", "1C_ideal_AU",
                 "1C_distorted_GC", "1C_distorted_AU"]),
    ("figS1.py", ["S1B", "S1C", "S1D"]),
    ("fig5.py", [("5C", "orig"), ("5C", "redo"),
                 ("5D", "orig"), ("5D", "redo")]),
    ("fig6.py", ["6A", "6B", "6C"]),
]


def isolate():
    """Delete every object that the panel left disabled.

    The panel functions enable exactly what belongs in the shot, so 'still
    disabled' is a reliable stand-in for 'not part of this panel'.
    """
    enabled = set(cmd.get_names("objects", enabled_only=1))
    for obj in cmd.get_names("objects"):
        if obj not in enabled:
            cmd.delete(obj)


def save_panel(script, panel, tag):
    cmd.reinitialize()
    cmd.do("run " + os.path.join(HERE, script))
    if isinstance(panel, tuple):
        cmd.do('panel("%s", "%s")' % panel)
    else:
        cmd.do('panel("%s")' % panel)
    isolate()
    path = os.path.join(OUT, tag + ".pse")
    cmd.save(path)
    n = len(cmd.get_names("objects"))
    print("  %-22s %d object(s)  %6.0f KB"
          % (tag + ".pse", n, os.path.getsize(path) / 1024.0))


def main():
    os.makedirs(OUT, exist_ok=True)
    for old in glob.glob(os.path.join(OUT, "*.pse")):
        os.remove(old)
    print("writing panels to %s\n" % OUT)
    for script, panels in JOBS:
        print(script)
        for p in panels:
            tag = "%s_%s" % p if isinstance(p, tuple) else p
            try:
                save_panel(script, p, tag)
            except Exception as exc:                       # noqa: BLE001
                print("  !! %s failed: %s" % (tag, exc))
    files = sorted(glob.glob(os.path.join(OUT, "*.pse")))
    print("\n%d .pse files in %s" % (len(files), OUT))


main()
