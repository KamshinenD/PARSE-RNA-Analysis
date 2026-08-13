"""PARSE -> PyMOL: highlight problematic base pairs during refinement.

    VENDORED COPY — the original lives in the engine repo at
    integrations/pymol/parse_pymol.py. It is copied here so Figure 6 can be
    regenerated from this repository alone. Re-copy it if the engine changes.

    Figure 6 is produced entirely by this plugin (see FIGURE6.md):
      panel A  parse_score 1GID, all
      panel B  parse_ideal on ; parse_goto A-235
      panel C  parse_ideal on ; parse_goto A-170

Load in PyMOL:    run /path/to/parse_pymol.py
Then:             parse_score               # score current object; worklist = Review tier
                  parse_score 6DN1, acceptable  # the minor-issue middle (polish pass)
                  parse_score 6DN1, all     # every non-Preferred pair
                  parse_clear               # remove the highlights

What it does each run:
  1. saves the object's nucleotides to a temp PDB (chains renamed to unique
     single characters so big multi-chain structures -- e.g. a ribosome with
     chain "1A" -- survive the legacy-PDB round-trip; res_ids are mapped back
     to the real chain names afterwards),
  2. runs the `parse` binary on it (sub-second, even on a ribosome),
  3. colors the whole structure gray (so unflagged pairs recede), then builds
     an overlay object `parse_overlay` of just the flagged base pairs, shown as
     sticks on a continuous yellow->red heat-map by Cerny severity (mild
     Allowed -> yellow, Of-Concern -> red),
  4. prints a ranked summary (incl. H-bond-count issues, which are summary-only).

Which pairs are flagged is gated by the whole-pair SCORE tier (`tier` arg to
parse_score): 'review' (score<75, default), 'acceptable' (75<=score<100), or
'all'. The score-tier decides *membership*; the Cerny |Z'| severity below still
decides how each flagged pair is *colored and ranked*.

Re-run after an edit/refinement round and the overlay rebuilds from scratch, so a
pair that's been fixed simply disappears. Only colors change — coordinates are
never modified (pass gray_background=0 to keep the object's existing colors).

Binary: set $PARSE_BINARY, or `parse_set_binary /path/to/parse`.
Severity coloring follows Cerny et al. (NAR 2026): Of Concern |Z'|>5, Allowed
|Z'|<=5. H-bond-count issues carry no ProSco/Z' tier and are reported in the
summary only, not colored.
"""
import json
import os
import re
import shutil
import string
import subprocess
import tempfile

# ---------------------------------------------------------------------------
# Pure logic (importable / testable without PyMOL)
# ---------------------------------------------------------------------------

#: issues that are not ProSco/Z'-scored -> excluded from coloring (summary only)
SUMMARY_ONLY = {"incorrect_hbond_count"}


def parse_res_id(res_id):
    """'A-5MU-54' -> ('A', '5MU', '54'); handles modified names and negative seq.

    Format is chain-resname-resseq. Only resseq may carry a sign, and resname
    never contains '-', so the first '-' splits chain, and the first '-' of the
    remainder splits resname from resseq (which keeps any leading '-').
    """
    i1 = res_id.index("-")
    chain = res_id[:i1]
    rest = res_id[i1 + 1:]
    i2 = rest.index("-")
    return chain, rest[:i2], rest[i2 + 1:]


def residue_selection(res_id):
    """PyMOL selection string for one PARSE res_id (negative resi escaped)."""
    chain, _name, seq = parse_res_id(res_id)
    resi = ("\\" + seq) if seq.startswith("-") else seq
    chain_tok = f"chain {chain} and " if chain else ""
    return f"({chain_tok}resi {resi})"


#: base ring + functional atoms — used to hide the cartoon "ladder" only where the
#: base is already drawn another way (focus sticks / context lines).
_BASE_ATOMS_SEL = "name N1+C2+N3+C4+C5+C6+N7+C8+N9+N6+O6+N2+O2+N4+O4"


def colored_issues(issue_details):
    """ProSco/Z'-scored issues only (drops summary-only ones like hbond count)."""
    return [d for d in (issue_details or []) if d.get("issue") not in SUMMARY_ONLY]


_ZPRIME_THRESHOLD = 5.0    # Cerny: |Z'| = 5 is the Allowed/Of-Concern boundary
_PROSCO_PREFERRED = 5.0    # Cerny: ProSco >= 5 is Preferred (not flagged)

#: Visualization-only down-weighting of issues that are NOT visibly obvious.
#: `hbond_angles` is a bond-geometry (chemistry) defect, not a shape defect -- a
#: pair can sit perfectly flat and well-twisted yet have strained H-bond angles,
#: so a high hbond_angles |Z'| would otherwise rocket a normal-looking pair to the
#: top of the worklist (and color it bright red). Multiplying its severity by < 1
#: keeps it flagged and reported with its TRUE sigma, but lets the visibly-
#: distorted pairs (propeller, opening, plane angle, shear, ...) dominate the
#: triage order. Affects ONLY the PyMOL worklist rank + highlight color -- not the
#: 0-100 score, the scorer, or any benchmark. (H-bond *distance* keeps full weight:
#: a long/short bond shows up as a visibly splayed pair.)
VISUAL_WEIGHT = {"hbond_angles": 0.3}


def rank_severity(detail):
    """Unclamped, visually-weighted Cerny severity of one parameter, for ranking.

    Base severity is |Z'| / 5 -- mirrors the scorer's `min(1, |Z'|/5)` (when
    ProSco < 5) but WITHOUT the cap, so pairs don't all saturate to 1.0 (a
    6.6-sigma deviation outranks a 5.1-sigma one). ProSco gates whether the
    parameter counts at all; |Z'| is the continuous magnitude. The result is then
    scaled by VISUAL_WEIGHT so non-shape issues (hbond_angles) rank/color lower.
    Returns 0 for Preferred parameters (ProSco >= 5).
    """
    ps, zp = detail.get("prosco"), detail.get("zprime")
    if ps is None or ps >= _PROSCO_PREFERRED:
        return 0.0
    w = VISUAL_WEIGHT.get(detail.get("issue"), 1.0)
    if zp is None:
        return w     # outside the distribution, no Z' magnitude available
    return w * abs(zp) / _ZPRIME_THRESHOLD


def worst_parameter(issue_details):
    """The colored issue with the largest unclamped |Z'| severity, or None."""
    ci = colored_issues(issue_details)
    if not ci:
        return None
    return max(ci, key=rank_severity)


def severity_to_rgb(sev):
    """Continuous yellow(0)->red(1) heat-map; clamps to [0,1]."""
    s = 0.0 if sev is None else max(0.0, min(1.0, float(sev)))
    return [1.0, 1.0 - s, 0.0]


# ---------------------------------------------------------------------------
# Pair-score tiers — whole-pair triage that GATES the worklist. Distinct from
# the per-parameter Cerny tiers (Preferred/Allowed/Of Concern) that COLOR it:
# the pair-score tier decides *which pairs a refiner reviews*, the |Z'| severity
# decides *how each flagged pair is colored/ranked*. Preferred = score 100,
# Acceptable = [75,100), Review = <75 (mirrors scoring/result.py::pair_tier).
# ---------------------------------------------------------------------------

def pair_tier(score):
    """Whole-pair triage tier from the 0-100 score, or None if unscored."""
    if score is None:
        return None
    if score >= 100.0:
        return "Preferred"
    if score >= 75.0:
        return "Acceptable"
    return "Review"


#: worklist tier filter word -> the pair-score tiers that populate the worklist.
_TIER_SETS = {
    "review": {"Review"},                       # default: the genuinely-bad pairs
    "acceptable": {"Acceptable"},               # the minor-issue middle (polish pass)
    "all": {"Acceptable", "Review"},            # everything non-Preferred
}


def tiers_for(name):
    """Resolve a tier-filter word to its tier set.

    Unknown/typo values fall back to the default 'review' set, so a mistyped
    filter never silently empties the worklist.
    """
    return _TIER_SETS.get(str(name).strip().lower(), _TIER_SETS["review"])


def pair_tier_of(p):
    """Pair-score tier of a scored pair dict (emitted `tier`, else derived)."""
    return p.get("tier") or pair_tier(p.get("score"))


def tier_counts(data):
    """Structure-level {Preferred, Acceptable, Review} counts.

    Prefers the binary's emitted `tier_summary`; falls back to deriving from
    per-pair scores for older `parse` outputs that predate the field.
    """
    ts = data.get("tier_summary")
    if ts:
        return ts
    counts = {"Preferred": 0, "Acceptable": 0, "Review": 0}
    for p in data.get("pairs", []):
        t = pair_tier_of(p)
        if t in counts:
            counts[t] += 1
    return counts


def pair_label(f):
    """Compact, readable name for a flagged pair, e.g. 'A226.U249 (chain A)'.

    The raw res_ids ('A-A-235' / 'A-U-239') run together when joined, so build
    'resnameSeq.resnameSeq' and append the chain (or both chains if they differ).
    """
    c1, n1, s1 = parse_res_id(f["res_id1"])
    c2, n2, s2 = parse_res_id(f["res_id2"])
    if c1 == c2:
        return f"{n1}{s1}·{n2}{s2} (chain {c1})"
    return f"{c1}/{n1}{s1}·{c2}/{n2}{s2}"


#: distinct color for H-bond-count pairs. Cool = a bonding/topology defect (a
#: missing or extra canonical H-bond), deliberately OFF the warm yellow->red
#: geometric-severity ramp so it never implies a shape distortion that isn't
#: there. Also far from the green ideal-overlay / in-range colors.
_HBOND_RGB = [0.15, 0.45, 1.0]


def _hbond_count_detail(p):
    """The pair's incorrect_hbond_count issue detail dict, or None."""
    for d in (p.get("issue_details") or []):
        if d.get("issue") == "incorrect_hbond_count":
            return d
    return None


def severity_headline(f):
    """One-line 'what's worst' string for a flagged pair.

    Geometry pairs: 'propeller 6.2σ' (worst parameter + |Z'|). H-bond-count
    pairs carry no |Z'|, so they headline as 'H-bond count' instead.
    """
    if f.get("kind") == "hbond":
        return "H-bond count"
    issue = f.get("worst_issue") or "?"
    z = f.get("worst_z")
    return f"{issue} {abs(z):.1f}σ" if z is not None else f"{issue} (no Z')"


def pair_score_str(f):
    """The pair's 0-100 PARSE quality score, or '?' when unavailable."""
    s = f.get("score")
    return f"{s:.0f}" if s is not None else "?"


def _pair_flag(p):
    """Classify one scored pair into a flag dict: geometry / hbond / clean.

    * ``geometry`` — has a ProSco/Z'-colorable parameter (shear, propeller,
      opening, ...); ranked/colored by its worst |Z'| on the warm ramp.
    * ``hbond``    — no colorable geometry, but a missing/extra canonical H-bond
      (incorrect_hbond_count); shown blue, ordered by its graded severity.
    * ``clean``    — nothing out of range (Preferred parameters only).

    Geometry wins when both are present (the count still shows in the breakdown).
    `rank` sorts geometry worst-first; `sev` (0..1) drives the color intensity.
    """
    base = {
        "res_id1": p["res_id1"], "res_id2": p["res_id2"],
        "sequence": p.get("sequence"),
        "lw_class": (p.get("classification") or {}).get("lw_class"),
        "pair_tier": pair_tier_of(p), "score": p.get("score"),
        "details": p.get("issue_details", []),
    }
    w = worst_parameter(p.get("issue_details"))
    if w is not None:
        rank = rank_severity(w)
        if rank > 0.0:
            base.update({
                "kind": "geometry", "rank": rank, "sev": min(1.0, rank),
                "worst_issue": w.get("issue"), "worst_z": w.get("zprime"),
                "tier": w.get("tier", ""), "issue": w.get("issue"),
            })
            return base
    hb = _hbond_count_detail(p)
    if hb is not None:
        sev = min(1.0, max(0.0, float(hb.get("severity") or 0.0)))
        base.update({
            "kind": "hbond", "rank": 0.0, "sev": sev,
            "worst_issue": "incorrect_hbond_count", "worst_z": None,
            "tier": "", "issue": "incorrect_hbond_count",
        })
        return base
    base.update({
        "kind": "clean", "rank": 0.0, "sev": 0.0,
        "worst_issue": None, "worst_z": None, "tier": "Preferred", "issue": None,
    })
    return base


def flagged_pairs(data, min_severity=0.0, tiers=None):
    """Pairs in the requested score-tier(s) with a colorable issue, worst-first.

    `tiers` is the set of pair-score tiers allowed into the worklist (default
    {"Review"}). A pair enters only if its pair-score tier is in `tiers` AND it
    has a colorable issue — either a ProSco/Z'-scored geometry parameter (warm
    ramp) OR a missing/extra canonical H-bond (blue). Geometry pairs are ranked
    worst-first by |Z'|; H-bond-count pairs (no |Z'|) are appended AFTER them,
    ordered by their graded severity. Returns _pair_flag dicts (see that fn).
    """
    if tiers is None:
        tiers = {"Review"}
    geom, hbond = [], []
    for p in data.get("pairs", []):
        if pair_tier_of(p) not in tiers:        # aggregate score-tier gate
            continue
        f = _pair_flag(p)
        if f["sev"] < min_severity:             # power-user color-intensity trim
            continue
        if f["kind"] == "geometry":
            geom.append(f)
        elif f["kind"] == "hbond":
            hbond.append(f)
        # clean -> not flagged
    geom.sort(key=lambda r: -r["rank"])
    hbond.sort(key=lambda r: -r["sev"])
    return geom + hbond


def hbond_count_pairs(data):
    """Pairs whose only/also-present issue is incorrect_hbond_count (summary)."""
    out = []
    for p in data.get("pairs", []):
        hb = [d for d in (p.get("issue_details") or [])
              if d.get("issue") == "incorrect_hbond_count"]
        if hb:
            out.append({"res_id1": p["res_id1"], "res_id2": p["res_id2"],
                        "sequence": p.get("sequence"),
                        "lw_class": (p.get("classification") or {}).get("lw_class")})
    return out


def format_summary(data, flagged, hbcount, obj, tier_label="review"):
    """Human-readable ranked summary string (ranked by ProSco/|Z'|, no score).

    Header shows the whole-structure pair-score tier counts, then the worklist
    line notes which tier filter it reflects (`tier_label`).
    """
    n_oc = sum(1 for f in flagged if f["tier"] == "Of Concern")
    n_al = sum(1 for f in flagged if f["tier"] == "Allowed")
    n_hb = sum(1 for f in flagged if f.get("kind") == "hbond")
    tc = tier_counts(data)
    lines = [f"PARSE: {obj}  {len(data.get('pairs', []))} pairs  "
             f"(structure: {tc.get('Preferred', 0)} Preferred · "
             f"{tc.get('Acceptable', 0)} Acceptable · {tc.get('Review', 0)} Review)",
             f"worklist [{tier_label} tier]: {len(flagged)} flagged "
             f"({n_oc} Of Concern, {n_al} Allowed, {n_hb} H-bond count)"]
    if flagged:
        lines.append("worst pairs (geometry by |Z'|, then H-bond-count in blue):")
        for f in flagged[:25]:
            extra = failing_params(f["details"]) or ("missing/extra canonical H-bond"
                                                      if f.get("kind") == "hbond" else "")
            lines.append(f"  {pair_label(f)}  "
                         f"{f.get('lw_class') or '?'} {f.get('sequence') or ''}  "
                         f"{severity_headline(f)}  [{extra}]")
    # Any remaining H-bond-count pairs outside the current tier (so not in the
    # worklist) — surfaced as text so they're never silently dropped.
    flagged_keys = {(f["res_id1"], f["res_id2"]) for f in flagged}
    only_count = [h for h in hbcount if (h["res_id1"], h["res_id2"]) not in flagged_keys]
    if only_count:
        lines.append(f"H-bond-count pairs outside the {tier_label} tier "
                     f"({len(only_count)}, not in this worklist):")
        for h in only_count[:25]:
            lines.append(f"  {pair_label(h)}  "
                         f"{h.get('lw_class') or '?'} {h.get('sequence') or ''}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# PyMOL integration
# ---------------------------------------------------------------------------

_BINARY = os.environ.get("PARSE_BINARY", "parse")
_OVERLAY = "parse_overlay"
_PROBLEMS = "parse_problems"
_FOCUS = "parse_focus"        # worklist: the current pair (bright sticks)
_CONTEXT = "parse_context"    # worklist: residues within a few A (faded lines)
_HB = "parse_hb"              # worklist: polar contacts within the current pair
_IDEAL = "parse_ideal_pair"   # ghost of the idealized pair, best-fit onto the focus
_IDEAL_HB = "parse_ideal_hb"  # the idealized pair's optimized H-bonds
_TMP_OBJ = "_parse_tmp"

#: worklist state set by parse_score, walked by parse_next / parse_prev / parse_goto
_QUEUE = []          # flagged pairs, ordered chain-grouped then worst-first
_POS = -1            # index of the current pair in _QUEUE (-1 = not started)
_SCORED_OBJ = None   # object the queue was built against
_STRUCT_SCORE = None # whole-structure scores from the last parse_score
_SCORED_DATA = None  # full scored JSON from the last parse_score (for parse_dump)
_SHOW_IDEAL = False  # sticky toggle: overlay the idealized pair on each focus
_LAST_FOCUS = None   # last pair-view rendered (so parse_ideal can act on it)
_CONTEXT_RADIUS = 8.0


def _plugin_dir():
    """Directory of this plugin file.

    PyMOL `run`s the script inside the `pymol` module namespace, so globals()
    __file__ points at PyMOL itself; PyMOL instead sets __script__ to this file.
    A plain import sets a correct __file__. Prefer whichever actually names this
    plugin; fall back to the cwd.
    """
    for var in ("__script__", "__file__"):
        p = globals().get(var)
        if p and "parse_pymol" in os.path.basename(p):
            return os.path.dirname(os.path.abspath(p))
    return os.getcwd()


#: idealized base-pair templates: <_IDEALS_DIR>/<lw_class>/<seq>.pdb. Defaults to
#: the repo's resources dir next to this plugin; override with PARSE_IDEALS_DIR or
#: the parse_set_ideals command.
_IDEALS_DIR = os.path.normpath(os.environ.get("PARSE_IDEALS_DIR") or os.path.join(
    _plugin_dir(), "..", "..", "resources", "basepair-idealized"))
#: single-char labels used to rename nucleotide chains before a legacy-PDB save
#: (A-Z, a-z, 0-9 = 62). Legacy PDB has a 1-char chain column and a 99999-atom
#: serial limit, so a big multi-char-chain structure (e.g. a ribosome with chain
#: "1A") would otherwise be saved with chains truncated/collapsed -> the scored
#: res_ids no longer match the in-PyMOL object and nothing gets highlighted.
_CHAIN_POOL = string.ascii_uppercase + string.ascii_lowercase + string.digits


def _resolve_binary():
    if os.path.sep in _BINARY and os.path.exists(_BINARY):
        return _BINARY
    found = shutil.which(_BINARY)
    if found:
        return found
    raise RuntimeError(
        f"parse binary not found ('{_BINARY}'). "
        f"Set it with: parse_set_binary /path/to/parse")


def _run_parse(pdb_path):
    binary = _resolve_binary()
    out = subprocess.run([binary, pdb_path, "--no-download"],
                         capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(f"parse failed:\n{out.stderr.strip()}")
    return json.loads(out.stdout)


def _backmap_chain(res_id, inv):
    """Rewrite a res_id's chain via the inverse remap (single-char -> real)."""
    chain, name, seq = parse_res_id(res_id)
    return f"{inv.get(chain, chain)}-{name}-{seq}"


def _score_object(cmd, obj):
    """Save `obj` to a temp PDB, score it, return the parsed JSON.

    The nucleotides are copied to a scratch object whose chains are renamed to
    unique single characters first, so the legacy-PDB round-trip can't truncate
    or collapse multi-character chains (e.g. a ribosome's "1A"/"1B"). The scored
    res_ids are then mapped back to the object's real chain names so the overlay
    selection matches the untouched original object. Protein/solvent is dropped
    (the scorer only scores base pairs); coordinates are never modified.

    Falls back to a plain whole-object save when there are no nucleotide chains
    or too many to fit the single-char pool.
    """
    # Broad nucleotide selection: PyMOL's nucleic polymer plus any residue that
    # carries a C1' sugar atom, so modified / HETATM bases are never dropped.
    nsel = f"(({obj}) and (polymer.nucleic or byres (({obj}) and name C1')))"
    chains = cmd.get_chains(nsel)

    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    tmp.close()
    inv = None
    try:
        if chains and len(chains) <= len(_CHAIN_POOL):
            fwd = {c: _CHAIN_POOL[i] for i, c in enumerate(chains)}
            inv = {v: k for k, v in fwd.items()}
            cmd.delete(_TMP_OBJ)
            cmd.create(_TMP_OBJ, nsel)
            try:
                for orig, single in fwd.items():
                    cmd.alter(f'{_TMP_OBJ} and chain "{orig}"', f"chain='{single}'")
                cmd.sort(_TMP_OBJ)
                cmd.save(tmp.name, _TMP_OBJ)
            finally:
                cmd.delete(_TMP_OBJ)
        else:
            # No nucleotides, or more chains than the single-char pool: best-effort
            # direct save (small single-char-chain structures already round-trip).
            cmd.save(tmp.name, obj)
        data = _run_parse(tmp.name)
    finally:
        os.unlink(tmp.name)

    if inv is not None:
        for p in data.get("pairs", []):
            p["res_id1"] = _backmap_chain(p["res_id1"], inv)
            p["res_id2"] = _backmap_chain(p["res_id2"], inv)
    data["pdb_id"] = obj   # the object's real name, not the scratch temp file
    return data


def parse_set_binary(path):
    """Point the plugin at a specific parse binary."""
    global _BINARY
    _BINARY = path
    print(f"PARSE: binary = {path}")


def parse_clear(_self=None):
    """Remove every PARSE overlay/selection; the original object is untouched."""
    from pymol import cmd
    global _QUEUE, _POS, _SCORED_OBJ
    if _SCORED_OBJ is not None:
        cmd.set("cartoon_transparency", 0.0, _SCORED_OBJ)
        cmd.show("cartoon", f"({_SCORED_OBJ}) and ({_BASE_ATOMS_SEL})")  # restore ladder
        _set_hud(cmd, "")
    for o in (_OVERLAY, _PROBLEMS, _FOCUS, _CONTEXT, _HB, _IDEAL, _IDEAL_HB):
        cmd.delete(o)
    _QUEUE, _POS, _SCORED_OBJ = [], -1, None
    print("PARSE: cleared")


def parse_score(obj=None, tier="review", min_severity=0.0, zoom=0,
                gray_background=1, _self=None):
    """Score the current (or named) object and highlight problematic base pairs.

    obj             object to score (default: first enabled object)
    tier            which pair-score tier(s) populate the worklist:
                      'review'     (default) score < 75  — the genuinely-bad pairs
                      'acceptable' 75 <= score < 100     — the minor-issue middle
                      'all'        both (everything non-Preferred)
                    Preferred (score 100) pairs are never flagged. The tier gates
                    *membership*; |Z'| severity still colors/ranks within it.
    min_severity    also drop flagged issues with severity < this (0..1; power-user
                    trim on top of the tier gate; default 0 = keep all in the tier)
    zoom            1 to zoom onto the flagged pairs
    gray_background 1 to color the whole structure gray (so unflagged geometry is
                    gray and only the flagged pairs carry color); 0 to leave the
                    object's existing colors. Coordinates are never modified.

    Examples:
        parse_score                 # current object, worklist = Review tier
        parse_score 6DN1, acceptable
        parse_score 6DN1, all, 0.3  # non-Preferred, min color severity 0.3
    """
    from pymol import cmd
    min_severity = float(min_severity)
    tier_label = str(tier).strip().lower()
    tiers = tiers_for(tier_label)

    if obj is None:
        objs = cmd.get_object_list()
        if not objs:
            print("PARSE: no object loaded")
            return
        obj = objs[0]

    try:
        data = _score_object(cmd, obj)
    except Exception as e:        # noqa: BLE001 - surface any failure to the user
        print(f"PARSE: error: {e}")
        return

    flagged = flagged_pairs(data, min_severity, tiers)
    hbcount = hbond_count_pairs(data)
    print(format_summary(data, flagged, hbcount, obj, tier_label))

    # Reset any previous worklist focus, then build a fresh queue for this run.
    global _QUEUE, _POS, _SCORED_OBJ, _STRUCT_SCORE, _SCORED_DATA
    cmd.delete(_FOCUS)
    cmd.delete(_CONTEXT)
    cmd.delete(_HB)
    cmd.set("cartoon_transparency", 0.0, obj)
    _QUEUE = build_queue(flagged)
    _POS = -1
    _SCORED_OBJ = obj
    _SCORED_DATA = data                    # full JSON, written out by parse_dump
    _STRUCT_SCORE = {                       # per-axis scores, shown by parse_overview
        "pairs":    data.get("pairs_score"),
        "residues": data.get("backbone_suiteness", data.get("backbone_score")),
        "n_pairs":  data.get("n_pairs"),
        "tiers":    tier_counts(data),
        "tier_label": tier_label,
    }

    # Gray the whole structure so Preferred geometry recedes and the flagged
    # pairs (colored sticks, built below) stand out. Color only, not coordinates.
    if int(gray_background):
        cmd.color("gray70", obj)

    # Rebuild a fresh overlay of just the flagged pairs.
    cmd.delete(_OVERLAY)
    cmd.delete(_PROBLEMS)
    if not flagged:
        tc = tier_counts(data)
        _set_hud(cmd, f"PARSE [{tier_label} tier]  0 flagged   "
                      f"(structure: {tc.get('Preferred', 0)} Preferred · "
                      f"{tc.get('Acceptable', 0)} Acceptable · {tc.get('Review', 0)} Review)")
        return
    _build_overlay(cmd, obj, flagged)

    cmd.select(_PROBLEMS, _OVERLAY)
    cmd.deselect()
    if int(zoom):
        cmd.zoom(_OVERLAY)

    # Guided HUD: bind nav keys (once) and show a starting status line.
    _bind_nav(cmd)
    n_oc = sum(1 for x in flagged if x["tier"] == "Of Concern")
    _set_hud(cmd, f"PARSE [{tier_label} tier]  {len(flagged)} flagged "
                  f"({n_oc} Of Concern)   press → to start the worklist   "
                  f"(Ctrl-I = inspect a clicked pair)")


def _build_overlay(cmd, obj, flagged):
    """Fast, non-destructive severity overlay.

    Instead of one giant residue selection (slow to parse) and a per-pair color
    loop (hundreds of PyMOL round-trips), we work on a single *copy*: stamp each
    flagged residue's clamped severity into the COPY's B-factor, trim the copy to
    the flagged residues, and color everything in one `spectrum` call. The B-factor
    is written only on the throwaway copy, so the original object's B-factors,
    colors and coordinates are never touched. ~150x faster on a ribosome.

    Two color classes: geometry residues carry b in [0,1] and get the warm
    yellow->red ramp; H-bond-count residues (bonding defect, no |Z'|) carry a
    negative sentinel and are colored a flat blue instead — geometry wins when a
    residue is in both. Non-flagged residues (b = -1) are removed.
    """
    # Geometry residues -> worst clamped severity (0..1). H-bond residues that are
    # NOT already geometry-flagged -> a keep-me sentinel (colored blue below).
    geo_sevmap, hb_res = {}, set()
    for f in flagged:
        for rid in (f["res_id1"], f["res_id2"]):
            ch, _, seq = parse_res_id(rid)
            if f.get("kind") == "hbond":
                hb_res.add((ch, seq))
            else:
                geo_sevmap[(ch, seq)] = max(geo_sevmap.get((ch, seq), 0.0), f["sev"])
    hb_res -= set(geo_sevmap)          # geometry wins the color for shared residues

    nsel = f"({obj}) and (polymer.nucleic or byres (({obj}) and name C1'))"
    cmd.create(_OVERLAY, nsel)
    cmd.alter(_OVERLAY,
              "b = geo.get((chain, resi), -0.25 if (chain, resi) in hb else -1.0)",
              space={"geo": geo_sevmap, "hb": hb_res})
    cmd.remove(f"{_OVERLAY} and b < -0.5")   # drop non-flagged; keep geo (>=0) + hb (-0.25)
    cmd.hide("everything", _OVERLAY)
    cmd.show("sticks", _OVERLAY)
    cmd.set("stick_radius", 0.18, _OVERLAY)
    # Geometry: clamped severity 0->1 maps yellow->red (matches severity_to_rgb).
    # (PyMOL's selection algebra has no ">="; the -0.1 cut splits geo>=0 from hb=-0.25.)
    cmd.spectrum("b", "yellow red", f"{_OVERLAY} and b > -0.1", minimum=0.0, maximum=1.0)
    # H-bond-count: flat blue, off the geometric ramp.
    cmd.set_color("parse_hbond_col", _HBOND_RGB)
    cmd.color("parse_hbond_col", f"{_OVERLAY} and b < -0.1")


# ---------------------------------------------------------------------------
# Worklist: walk the flagged pairs one at a time (parse_next / prev / goto)
# ---------------------------------------------------------------------------

def build_queue(flagged):
    """Order flagged pairs chain-grouped, worst-first within each chain.

    Geometry pairs are grouped by chain (first-seen order = worst-chain first),
    worst-first within each chain, so walking the queue stays in one region
    instead of teleporting. H-bond-count pairs (a subtler bonding defect, no
    |Z'|) are kept OUT of the chain grouping and appended at the very bottom,
    in their incoming severity order — so the worklist ends with them.
    """
    geom = [f for f in flagged if f.get("kind") != "hbond"]
    hbond = [f for f in flagged if f.get("kind") == "hbond"]
    by_chain, order = {}, []
    for f in geom:
        ch = parse_res_id(f["res_id1"])[0]
        if ch not in by_chain:
            by_chain[ch] = []
            order.append(ch)
        by_chain[ch].append(f)
    queue = []
    for ch in order:
        queue.extend(by_chain[ch])
    queue.extend(hbond)                 # H-bond-count pairs at the very bottom
    return queue


def failing_params(details, max_n=4):
    """Short 'which way is it distorted' string from the worst colored issues.

    e.g. 'shear Z'=+11.3, opening Z'=+10.3, propeller Z'=+8.5'. Of-Concern first,
    then by |Z'|; summary-only issues (hbond count) are excluded.
    """
    ci = [d for d in (details or []) if d.get("issue") not in SUMMARY_ONLY]
    ci.sort(key=lambda d: (d.get("tier") != "Of Concern",
                           -abs(d.get("zprime") or 0.0)))
    parts = []
    for d in ci[:max_n]:
        z = d.get("zprime")
        parts.append(f"{d['issue']} Z'={z:+.1f}" if z is not None else d["issue"])
    return ", ".join(parts)


def phenix_selection(f):
    """Phenix-style atom-selection string for a flagged pair (for hand-off)."""
    c1, _, s1 = parse_res_id(f["res_id1"])
    c2, _, s2 = parse_res_id(f["res_id2"])
    return f"(chain {c1} and resseq {s1}) or (chain {c2} and resseq {s2})"


def _focus(cmd, idx, dim=1):
    """Frame worklist queue[idx] and move the worklist cursor there."""
    global _POS
    _POS = idx
    _render_focus(cmd, _QUEUE[idx], dim=dim, pos=idx + 1, total=len(_QUEUE))


def _render_focus(cmd, f, dim=1, pos=None, total=None):
    """Frame one pair-view dict `f`: bright pair + faded context, H-bonds, label.

    `f` is a flagged_pairs-style dict (use _pair_view for a clean pair). With
    `pos`/`total` it's a worklist stop ('[3/45]'); without them it's an ad-hoc
    residue lookup (parse_goto A-169) and the worklist cursor is left untouched.
    A pair with nothing out of range is drawn green and reported as in-range.
    """
    global _LAST_FOCUS
    _LAST_FOCUS = f
    obj = _SCORED_OBJ
    flagged = f.get("worst_issue") is not None
    s1, s2 = residue_selection(f["res_id1"]), residue_selection(f["res_id2"])
    pair_sel = f"({obj}) and ({s1} or {s2})"

    cmd.delete(_FOCUS)
    cmd.delete(_CONTEXT)
    cmd.delete(_HB)

    # De-emphasise the global view so the current pair reads clearly.
    cmd.disable(_OVERLAY)
    if int(dim):
        cmd.set("cartoon_transparency", 0.85, obj)

    # The pair itself: bright sticks. Blue for an H-bond-count defect, warm
    # severity color for a geometry defect, green when nothing is out of range.
    cmd.create(_FOCUS, pair_sel)
    cmd.hide("everything", _FOCUS)
    cmd.show("sticks", _FOCUS)
    cmd.set("stick_radius", 0.25, _FOCUS)
    if f.get("kind") == "hbond":
        focus_col = _HBOND_RGB
    elif flagged:
        focus_col = severity_to_rgb(f["sev"])
    else:
        focus_col = [0.35, 0.78, 0.35]
    cmd.set_color("parse_focus_col", focus_col)
    cmd.color("parse_focus_col", _FOCUS)

    # Local context: nearby residues as thin gray lines (what's crowding it).
    cmd.create(_CONTEXT, f"byres (({obj}) within {_CONTEXT_RADIUS} of ({pair_sel}))"
                         f" and not ({s1} or {s2})")
    cmd.hide("everything", _CONTEXT)
    cmd.show("lines", _CONTEXT)
    cmd.color("gray50", _CONTEXT)

    # Drop the cartoon base "ladder" ONLY where the base is already drawn (focus
    # sticks + context lines) -- it's redundant there. Reset first (re-show all
    # base cartoon), then hide it for this focus+context region; residues outside
    # the context keep their ladder as their only base representation. The
    # backbone cartoon is untouched, so the trace stays continuous.
    cmd.show("cartoon", f"({obj}) and ({_BASE_ATOMS_SEL})")
    region = (f"({obj}) and ({_BASE_ATOMS_SEL}) and ({s1} or {s2} or "
              f"byres (({obj}) within {_CONTEXT_RADIUS} of ({pair_sel})))")
    cmd.hide("cartoon", region)

    # Polar contacts within the pair (a stand-in for its H-bonds).
    cmd.distance(_HB, _FOCUS, _FOCUS, mode=2)
    cmd.hide("labels", _HB)

    # One label on the pair: position tag, pair, class, score, and status.
    tag = f"#{pos} " if pos else ""
    headline = severity_headline(f) if flagged else "in range"
    label = (f"{tag}{pair_label(f)} "
             f"{f.get('lw_class') or '?'} {f.get('sequence') or ''} "
             f"— score {pair_score_str(f)}/100 ({headline})")
    cmd.label(f"{_FOCUS} and {s1} and name C1'", repr(label))

    cmd.zoom(_FOCUS, 5)

    # Console: score + which parameters are out of range (or that none are).
    where = f"[{pos}/{total}]" if pos else "(lookup)"
    kind = f.get("kind")
    if kind == "hbond":
        tier, status = "H-bond count", "worst: H-bond count (bonding defect)"
    elif flagged:
        tier = f.get("tier") or ""
        status = f"worst: {severity_headline(f)}"
    else:
        tier, status = "Preferred", "all parameters in range"
    print(f"PARSE {where} {pair_label(f)}  "
          f"{f.get('lw_class') or '?'} {f.get('sequence') or ''}  "
          f"score {pair_score_str(f)}/100  {status}  ({tier})")
    if kind == "hbond":
        print("  out of range: missing or extra canonical H-bond (geometry in range)")
    elif flagged:
        print(f"  out of range: {failing_params(f['details']) or '(none colored)'}")
    else:
        print("  out of range: none — every parameter is within Preferred range")
    # Phenix selection hidden for now (re-enable when needed):
    # print(f"  phenix:     {phenix_selection(f)}")

    # On-screen HUD (viewport title): always-visible position + current pair.
    where_hud = f"{pos}/{total}" if pos else "lookup"
    if kind == "hbond":
        hud_status = "H-bond count (blue)"
    elif flagged:
        hud_status = f"{severity_headline(f)} [{tier}]"
    else:
        hud_status = f"in range [{tier}]"
    _set_hud(cmd, f"PARSE {where_hud}  {pair_label(f)}  "
                  f"{f.get('lw_class') or '?'} {f.get('sequence') or ''}  "
                  f"score {pair_score_str(f)}/100  {hud_status}   "
                  f"(→ next  ← prev  Ctrl-I = inspect click)")

    # Idealized-pair ghost (target geometry), only while the sticky toggle is on.
    if _SHOW_IDEAL:
        if not _overlay_ideal(cmd, f):
            print(f"  ideal: no template for {f.get('lw_class')} {f.get('sequence')}")
    else:
        _clear_ideal(cmd)


def _ideal_template_path(f):
    """Path to the idealized-pair PDB for pair-view `f`, or None if absent.

    Templates live at <_IDEALS_DIR>/<lw_class>/<seq>.pdb (e.g. cWW/GC.pdb). The
    first class of an ambiguous label is used; sequence is taken from the dump.
    """
    lw = (f.get("lw_class") or "").split("|")[0].strip()
    seq = (f.get("sequence") or "").strip()
    if not lw or not seq:
        return None
    # Try the sequence as-is, then reversed: symmetric classes (cWW, tWW, ...) are
    # stored under one canonical order, so a G·C pair must fall back to CG.pdb.
    for s in (seq, seq[::-1]):
        path = os.path.join(_IDEALS_DIR, lw, f"{s}.pdb")
        if os.path.exists(path):
            return path
    return None


def _clear_ideal(cmd):
    """Remove the idealized-pair ghost objects (idempotent)."""
    cmd.delete(_IDEAL)
    cmd.delete(_IDEAL_HB)


def _fit_pairs(cmd, ideal_res_sel, actual_res_sel):
    """Alternating (ideal, actual) per-atom selections for pair_fit.

    One pair per atom name present in BOTH residues, so pair_fit superposes
    least-squares over the shared base atoms (ideal is the mobile object).
    """
    ni = {a.name for a in cmd.get_model(ideal_res_sel).atom}
    na = {a.name for a in cmd.get_model(actual_res_sel).atom}
    pairs = []
    for nm in sorted(ni & na):
        pairs.append(f"({ideal_res_sel}) and name {nm}")
        pairs.append(f"({actual_res_sel}) and name {nm}")
    return pairs


def _overlay_ideal(cmd, f):
    """Best-fit the idealized pair for `f` onto the actual pair; draw it green.

    Returns True if drawn, False if there's no template or too few shared atoms
    for a stable fit. The fit is least-squares over every atom name shared by
    ideal & actual in BOTH residues, so both bases show their deviation.
    """
    _clear_ideal(cmd)
    path = _ideal_template_path(f)
    if path is None or _SCORED_OBJ is None:
        return False
    obj = _SCORED_OBJ
    cmd.load(path, _IDEAL, zoom=0)            # zoom=0: don't fly the camera to the
                                             # template's origin coords before we fit
    # Match each actual residue to the ideal chain of the SAME base by shared atom
    # count -- handles templates stored in the reverse sequence order (GC<->CG),
    # where ideal chain A is the partner of res_id1, not res_id1 itself.
    a1 = f"({obj}) and {residue_selection(f['res_id1'])}"
    a2 = f"({obj}) and {residue_selection(f['res_id2'])}"
    iA = f"({_IDEAL}) and chain A and resi 1"
    iB = f"({_IDEAL}) and chain B and resi 2"
    direct = _fit_pairs(cmd, iA, a1) + _fit_pairs(cmd, iB, a2)
    swapped = _fit_pairs(cmd, iB, a1) + _fit_pairs(cmd, iA, a2)
    pairs = direct if len(direct) >= len(swapped) else swapped
    if len(pairs) < 8:                       # < 4 shared atoms -> unreliable fit
        _clear_ideal(cmd)
        return False
    try:
        cmd.pair_fit(*pairs)
    except Exception:                        # noqa: BLE001 - bail on a fit failure
        _clear_ideal(cmd)
        return False
    # Green ghost: thin, semi-transparent sticks, distinct from the real pair.
    cmd.hide("everything", _IDEAL)
    cmd.show("sticks", _IDEAL)
    cmd.set("stick_radius", 0.12, _IDEAL)
    cmd.set("stick_transparency", 0.3, _IDEAL)
    cmd.color("green", _IDEAL)
    # The ideal's optimized H-bonds (green dashes).
    cmd.distance(_IDEAL_HB, _IDEAL, _IDEAL, mode=2)
    cmd.hide("labels", _IDEAL_HB)
    cmd.color("green", _IDEAL_HB)
    return True


def _set_hud(cmd, text):
    """Show a one-line status in the viewport (object title), or clear it."""
    if _SCORED_OBJ is None:
        return
    try:
        cmd.set_title(_SCORED_OBJ, 1, text)
    except Exception:        # noqa: BLE001 - title is cosmetic, never fail on it
        pass


def _require_queue():
    if not _QUEUE:
        print("PARSE: no worklist - run parse_score first")
        return False
    return True


def parse_next(_self=None):
    """Step to the next (less severe) flagged pair in the worklist."""
    from pymol import cmd
    if not _require_queue():
        return
    if _POS >= len(_QUEUE) - 1:
        print(f"PARSE: at the end of the worklist ({len(_QUEUE)} pairs)")
        return
    _focus(cmd, _POS + 1)


def parse_prev(_self=None):
    """Step back to the previous (more severe) flagged pair."""
    from pymol import cmd
    if not _require_queue():
        return
    if _POS <= 0:
        print("PARSE: at the start of the worklist")
        return
    _focus(cmd, _POS - 1)


def parse_overview(_self=None):
    """Zoom back out to the whole-structure overview after parse_next/parse_prev.

    Undims the cartoon, drops the per-pair focus objects, re-shows the overlay
    of all flagged pairs, resets the worklist to the top, and frames the whole
    structure. Does NOT re-score (use parse_score to recompute).
    """
    from pymol import cmd
    if _SCORED_OBJ is None:
        print("PARSE: nothing scored - run parse_score first")
        return
    global _POS
    for o in (_FOCUS, _CONTEXT, _HB, _IDEAL, _IDEAL_HB):   # drop per-pair objects
        cmd.delete(o)
    cmd.set("cartoon_transparency", 0.0, _SCORED_OBJ)   # undim the structure
    cmd.show("cartoon", f"({_SCORED_OBJ}) and ({_BASE_ATOMS_SEL})")  # restore base ladder
    cmd.enable(_OVERLAY)                        # re-show all flagged pairs
    _POS = -1                                   # overview = not on any one pair
    cmd.zoom(_SCORED_OBJ)                       # frame the whole structure

    # Whole-RNA score line (console + HUD).
    def _f(x):
        return f"{x:.1f}" if isinstance(x, (int, float)) else "?"
    def _f3(x):        # backbone_suiteness is a 0-1 fraction (wwPDB convention)
        return f"{x:.3f}" if isinstance(x, (int, float)) else "?"
    ss = _STRUCT_SCORE or {}
    tc = ss.get("tiers") or {}
    tlabel = ss.get("tier_label", "review")
    score_line = (f"{_SCORED_OBJ}: PARSE pairs {_f(ss.get('pairs'))}/100, "
                  f"backbone suiteness {_f3(ss.get('residues'))}  "
                  f"(structure: {tc.get('Preferred', 0)} Preferred · "
                  f"{tc.get('Acceptable', 0)} Acceptable · {tc.get('Review', 0)} Review)"
                  f"  — {len(_QUEUE)} flagged [{tlabel} tier]")
    print(f"PARSE: overview — {score_line}")
    _set_hud(cmd, f"PARSE overview   {score_line}   (→ parse_next to step through)")


def _parse_residue_token(token):
    """Parse a residue reference into (chain, resseq), or None.

    Accepts 'A-169', 'A/169', 'A 169', 'A:169', 'A169', or a full PARSE res-id
    'chain-resname-resseq' (e.g. 'A-G-169', '1A-5MU-54'). Chain may be multi-char.
    """
    t = token.strip()
    if t.count("-") >= 2:                      # looks like a full res-id
        try:
            chain, _name, seq = parse_res_id(t)
            return chain, seq
        except (ValueError, IndexError):
            pass
    m = re.match(r"^([A-Za-z0-9]+?)[-/ :]+(-?\d+)$", t)   # chain<sep>resseq
    if m:
        return m.group(1), m.group(2)
    m = re.match(r"^([A-Za-z]+)(-?\d+)$", t)              # chainresseq, no sep
    if m:
        return m.group(1), m.group(2)
    return None


def _pair_view(p):
    """Build a flag view dict from a raw scored pair `p` (for ad-hoc lookups).

    Delegates to _pair_flag, so a geometry pair is drawn on the warm ramp, an
    H-bond-count pair blue, and a clean (Preferred) pair green / 'in range'.
    """
    return _pair_flag(p)


def _find_scored_pairs(chain, resseq):
    """All scored pairs (flagged or clean) containing residue chain+resseq."""
    found = []
    for p in (_SCORED_DATA or {}).get("pairs", []):
        for rid in (p.get("res_id1"), p.get("res_id2")):
            if not rid:
                continue
            try:
                c, _n, s = parse_res_id(rid)
            except (ValueError, IndexError):
                continue
            if c == chain and s == resseq:
                found.append(p)
                break
    return found


def parse_goto(which, _self=None):
    """Jump to a pair by worklist rank, or inspect any residue by chain+number.

    - `parse_goto 12`     -> the 12th worklist entry (flagged, worst-first).
    - `parse_goto A-169`  -> the pair containing chain A residue 169, EVEN IF it
                             is a clean (Preferred) pair not in the worklist; prints
                             its score and which parameters (if any) are out of range.
    Residue forms: 'A-169', 'A/169', 'A 169', 'A169', or a full 'A-G-169'.
    """
    from pymol import cmd
    if not _require_queue():
        return
    token = str(which).strip()
    if token.isdigit():
        idx = int(token) - 1
        if not 0 <= idx < len(_QUEUE):
            print(f"PARSE: out of range 1..{len(_QUEUE)}")
            return
        _focus(cmd, idx)
        return

    res = _parse_residue_token(token)
    if res is None:
        print(f"PARSE: '{token}' is not a rank or a chain+number (e.g. 12 or A-169)")
        return
    chain, resseq = res

    # A flagged worklist pair -> normal focus (keeps its rank + advances cursor).
    for i, f in enumerate(_QUEUE):
        for rid in (f["res_id1"], f["res_id2"]):
            c, _n, s = parse_res_id(rid)
            if c == chain and s == resseq:
                _focus(cmd, i)
                return

    # Otherwise look across ALL scored pairs (clean ones included).
    hits = _find_scored_pairs(chain, resseq)
    if not hits:
        print(f"PARSE: {chain}/{resseq} is not in any scored pair "
              f"(unpaired, or not in this structure)")
        return
    _render_focus(cmd, _pair_view(hits[0]), dim=1)       # ad-hoc: cursor untouched
    if len(hits) > 1:
        partners = ", ".join(pair_label(_pair_view(p)) for p in hits[1:])
        print(f"  note: {chain}/{resseq} is also in {len(hits) - 1} other pair(s): {partners}")


def parse_list(n=25, _self=None):
    """Print the ranked worklist (rank, pair, class, worst |Z'|, distortion)."""
    if not _require_queue():
        return
    n = int(n)
    print(f"PARSE worklist: {len(_QUEUE)} pairs (geometry by |Z'|, then H-bond-count)")
    for i, f in enumerate(_QUEUE[:n]):
        mark = ">" if i == _POS else " "
        note = failing_params(f["details"], 2) or ("H-bond count"
                                                   if f.get("kind") == "hbond" else "")
        print(f" {mark}{i + 1:>4}  {pair_label(f)}  "
              f"{f.get('lw_class') or '?'} {f.get('sequence') or ''}  "
              f"score {pair_score_str(f)}/100  "
              f"{severity_headline(f)} {f['tier']}  [{note}]")
    if len(_QUEUE) > n:
        print(f"  ... {len(_QUEUE) - n} more (parse_list {len(_QUEUE)} for all)")


def parse_dump(path=None, _self=None):
    """Write the current scored data to a JSON file the user can inspect.

    This is exactly the JSON the last parse_score computed (every pair with its
    score, issues, geometry and rigid-body params, plus the structure-level
    scores) -- including any coordinate edits you re-scored. Nothing is written
    unless you call this. Default path: '<object>_pairs.json' in the current
    directory; pass any name or path to put it wherever you like.
    """
    import json
    import os
    if _SCORED_DATA is None or _SCORED_OBJ is None:
        print("PARSE: nothing scored - run parse_score first")
        return
    path = str(path) if path else f"{_SCORED_OBJ}_pairs.json"
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    try:
        with open(path, "w") as fh:
            json.dump(_SCORED_DATA, fh, indent=1)
    except OSError as e:        # noqa: BLE001 - report write failures to the user
        print(f"PARSE: could not write {path}: {e}")
        return
    n = len(_SCORED_DATA.get("pairs", []))
    print(f"PARSE: wrote {n} scored pair(s) to {path}")


def parse_load(path, tier="review", _self=None):
    """Rebuild the worklist from a parse_dump JSON -- WITHOUT re-scoring.

    Lets anyone (e.g. a colleague you sent the .pse + JSON to) step through a
    previously scored structure -- parse_next / parse_prev / parse_list /
    parse_overview / parse_dump -- with no parse binary. The object named
    in the JSON (its 'pdb_id') should already be loaded, e.g. from the .pse.
    `tier` selects which pair-score tier(s) populate the worklist (default
    'review'; see parse_score).
    """
    from pymol import cmd
    import json
    try:
        with open(str(path)) as fh:
            data = json.load(fh)
    except (OSError, ValueError) as e:        # noqa: BLE001 - report bad file
        print(f"PARSE: could not read {path}: {e}")
        return
    obj = data.get("pdb_id")
    if not obj:
        print(f"PARSE: {path} has no 'pdb_id' - not a parse_dump file?")
        return
    if obj not in cmd.get_object_list():
        print(f"PARSE: note - object '{obj}' is not loaded; navigation needs it "
              f"(open its .pse/structure first).")

    tier_label = str(tier).strip().lower()
    global _QUEUE, _POS, _SCORED_OBJ, _STRUCT_SCORE, _SCORED_DATA
    _SCORED_DATA = data
    _SCORED_OBJ = obj
    _QUEUE = build_queue(flagged_pairs(data, tiers=tiers_for(tier_label)))
    _POS = -1
    _STRUCT_SCORE = {
        "pairs":    data.get("pairs_score"),
        "residues": data.get("backbone_suiteness", data.get("backbone_score")),
        "n_pairs":  data.get("n_pairs"),
        "tiers":    tier_counts(data),
        "tier_label": tier_label,
    }
    _bind_nav(cmd)
    print(f"PARSE: loaded {len(_QUEUE)} flagged pair(s) [{tier_label} tier] for "
          f"'{obj}' from {path}  (parse_next to step through — no re-scoring)")


def parse_ideal(state="toggle", _self=None):
    """Overlay the idealized version of the current pair (green ghost) to show
    the target geometry — what the pair *should* look like, so a refiner can see
    how to fix it.

        parse_ideal on     keep the ideal overlaid on every parse_next/parse_goto
        parse_ideal off    stop overlaying it
        parse_ideal        toggle for the current pair

    The ideal is best-fit onto the actual pair (least-squares over the shared base
    atoms), so both bases show how far they sit from ideal. Templates come from
    resources/basepair-idealized/<class>/<seq>.pdb (set with parse_set_ideals).
    """
    from pymol import cmd
    global _SHOW_IDEAL
    s = str(state).lower()
    if s in ("off", "0", "false", "hide", "no"):
        _SHOW_IDEAL = False
    elif s in ("on", "1", "true", "show", "yes"):
        _SHOW_IDEAL = True
    else:
        _SHOW_IDEAL = not _SHOW_IDEAL        # bare `parse_ideal` toggles

    if not _SHOW_IDEAL:
        _clear_ideal(cmd)
        print("PARSE: ideal overlay OFF")
        return
    print("PARSE: ideal overlay ON — green ghost = ideal geometry for this pair")
    if _LAST_FOCUS is None:
        print("  (step to a pair with parse_next / parse_goto to see it)")
        return
    f = _LAST_FOCUS
    if not _overlay_ideal(cmd, f):
        print(f"  ideal: no template for {f.get('lw_class')} {f.get('sequence')}")
        return
    # Frame the ghost together with the real pair (covers the parse_overview case,
    # where no pair is focused so the camera would otherwise stay zoomed out).
    s1, s2 = residue_selection(f["res_id1"]), residue_selection(f["res_id2"])
    cmd.zoom(f"({_IDEAL}) or (({_SCORED_OBJ}) and ({s1} or {s2}))", 5)


def parse_set_ideals(path, _self=None):
    """Set the directory holding idealized base-pair templates (<class>/<seq>.pdb)."""
    global _IDEALS_DIR
    _IDEALS_DIR = os.path.normpath(str(path))
    ok = os.path.isdir(_IDEALS_DIR)
    print(f"PARSE: idealized-pair dir = {_IDEALS_DIR}" + ("" if ok else "  (NOT FOUND)"))


def parse_info(selection="pk1", _self=None):
    """Inspect a clicked/selected residue: jump the worklist to its pair.

    Click an atom in the viewport (that makes the `pk1` selection), then run
    `parse_info` (or press `i`). If the residue is a flagged pair, the worklist
    jumps there and prints its breakdown; otherwise it says so.
    """
    from pymol import cmd
    if not _require_queue():
        return
    picks = []
    try:
        cmd.iterate(selection, "picks.append((chain, resi))", space={"picks": picks})
    except Exception:        # noqa: BLE001 - empty/invalid selection
        picks = []
    if not picks:
        print("PARSE: nothing selected — click an atom first, or pass a selection")
        return
    # A selection can span several residues (a whole pair, a box-select, a named
    # `sele`); check every distinct residue in it, not just the first atom, and
    # jump to the first one that belongs to a flagged pair.
    wanted, seen = [], set()
    for chain, resi in picks:
        if (chain, resi) not in seen:
            seen.add((chain, resi))
            wanted.append((chain, resi))
    for i, f in enumerate(_QUEUE):
        for rid in (f["res_id1"], f["res_id2"]):
            c, _name, s = parse_res_id(rid)
            if (c, s) in seen:
                _focus(cmd, i)
                return
    shown = ", ".join(f"{c}/{r}" for c, r in wanted[:6])
    print(f"PARSE: {shown} — not in the worklist (Preferred or unscored). "
          f"`parse_list` shows the {len(_QUEUE)} flagged pair(s).")


#: nav keybindings — PyMOL's set_key only takes special/modifier keys (not plain
#: letters), so arrows drive the worklist and Ctrl-I inspects the clicked atom.
_NAV_KEYS = (("right", lambda *a: parse_next()),
             ("left", lambda *a: parse_prev()),
             ("CTRL-I", lambda *a: parse_info("pk1")))


def _bind_nav(cmd):
    """Bind → next / ← prev / Ctrl-I inspect (idempotent, never raises)."""
    for key, fn in _NAV_KEYS:
        try:
            cmd.set_key(key, fn)
        except Exception:        # noqa: BLE001 - some keys unmappable on a platform
            pass


def parse_keys(state="on", _self=None):
    """Bind/unbind worklist keys: → next, ← prev, Ctrl-I inspect clicked atom."""
    from pymol import cmd
    if str(state).lower() in ("off", "0", "false"):
        for key, _fn in _NAV_KEYS:
            try:
                cmd.set_key(key, lambda *a: None)
            except Exception:    # noqa: BLE001
                pass
        print("PARSE: keyboard navigation off")
        return
    _bind_nav(cmd)
    print("PARSE keys: → next   ← prev   Ctrl-I = inspect clicked atom "
          "(or click + `parse_info`).  parse_keys off to unbind.")


# Register as PyMOL commands (no-op if imported outside PyMOL).
try:
    from pymol import cmd as _cmd
    _cmd.extend("parse_score", parse_score)
    _cmd.extend("parse_clear", parse_clear)
    _cmd.extend("parse_set_binary", parse_set_binary)
    _cmd.extend("parse_next", parse_next)
    _cmd.extend("parse_prev", parse_prev)
    _cmd.extend("parse_overview", parse_overview)
    _cmd.extend("parse_goto", parse_goto)
    _cmd.extend("parse_list", parse_list)
    _cmd.extend("parse_dump", parse_dump)
    _cmd.extend("parse_load", parse_load)
    _cmd.extend("parse_info", parse_info)
    _cmd.extend("parse_keys", parse_keys)
    _cmd.extend("parse_ideal", parse_ideal)
    _cmd.extend("parse_set_ideals", parse_set_ideals)
    print("PARSE plugin loaded: parse_score / parse_clear / parse_set_binary / "
          "parse_next / parse_prev / parse_overview / parse_goto / parse_list / "
          "parse_dump / parse_load / parse_info / parse_keys")
except ImportError:
    pass
