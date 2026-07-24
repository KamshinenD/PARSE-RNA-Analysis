"""Canonical orientation for self-reciprocal LW classes.

A pair like G-C/cWW and the reverse C-G/cWW are the same physical pair
viewed from either end — for self-reciprocal LW classes (where both
edges are the same kind), the only difference is which residue is
listed first. Without this module, the build pipeline splits such
pairs into two mirror cells, halving the data available to each.

Canonical form: residues in alphabetical order. So G-C/cWW becomes
C-G/cWW, U-A/cWW becomes A-U/cWW, etc. Heteroreciprocal LW classes
(cWH, cWS, tHS, …) are NOT canonicalized — their two orderings are
genuinely different pair types.

Some pair-geometry parameters flip sign when the residues are swapped
(shear, stagger, buckle) — the script that builds distributions and the
runtime that queries them both apply this flip before the |Z'| math.

Used by:
  - scripts/reference/build_distributions.py
  - scripts/reference/calculate_penalty_weights.py
  - src/pair_finder/scoring/issues.py (runtime)
"""

from __future__ import annotations

# LW classes where the two edges are the same kind (W·W, H·H, S·S).
# These collapse to one canonical orientation.
SELF_RECIPROCAL_LW: frozenset[str] = frozenset({
    "cWW", "tWW", "cHH", "tHH", "cSS", "tSS",
})

# Pair-geometry parameters whose sign flips under residue swap.
# Tsukuba convention: shear (lateral), stagger (out-of-plane), buckle
# (rotation about long axis) are all directional. stretch, propeller,
# opening are reflection-symmetric and stay the same on swap.
SIGN_FLIPPING_PAIR_PARAMS: frozenset[str] = frozenset({
    "shear", "stagger", "buckle",
})


def canonicalize_bp(bp_type: str | None, lw_class: str | None) -> tuple[str | None, bool]:
    """Return (canonical_bp_type, swapped).

    For self-reciprocal LW classes, canonicalize by alphabetical base order:
    G-C/cWW → C-G/cWW (swapped=True), U-A/cWW → A-U/cWW (swapped=True).
    For non-self-reciprocal classes or non-hyphen bp_types ('_ANY'),
    returns the input unchanged with swapped=False.
    """
    if bp_type is None or lw_class is None:
        return bp_type, False
    if lw_class not in SELF_RECIPROCAL_LW:
        return bp_type, False
    if "-" not in bp_type:
        return bp_type, False
    left, right = bp_type.split("-", 1)
    if left <= right:
        return bp_type, False
    return f"{right}-{left}", True


def maybe_flip(value: float | None, source: str, param: str, swapped: bool) -> float | None:
    """If canonicalization required a residue swap AND this is an
    orientation-dependent pair-geometry parameter, flip the sign.

    H-bond parameters and stretch/propeller/opening are unaffected.
    """
    if value is None or not swapped:
        return value
    if source == "pairs" and param in SIGN_FLIPPING_PAIR_PARAMS:
        return -value
    return value
