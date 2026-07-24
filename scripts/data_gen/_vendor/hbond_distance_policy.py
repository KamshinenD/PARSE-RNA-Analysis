"""Policy for which H-bond distances define a base pair's distance score.

A pair must carry >= 1 STRONG base-base bond (heavy-atom donor-acceptor
distance <= MAX_STRONG_DISTANCE) to be accepted by the finder (except SS/HH,
which may be all-weak via the geometry/template path). This module is the
single source of truth for which bond distances feed the distance
distribution (build side) and which are scored at runtime:

  - Watson-Crick cells (C-G/cWW, A-U/cWW): ALL base-base bonds. WC geometry is
    rigid, so a long bond anywhere in the pair is anomalous and penalized.
  - Every other cell: only the STRONG (primary) bond = the strongest
    (min-distance) base-base bond, and only if it is <= MAX_STRONG_DISTANCE.
    Weaker secondary bonds are exempt from the distance penalty (they still
    count toward num_base_hbonds via the incorrect_hbond_count parameter).
    A pair with no strong bond (SS/HH geometry-exception acceptance)
    contributes no distance and pays no distance penalty.

Sub-physical contacts (< MIN_PHYSICAL_DISTANCE, e.g. ~1.3 A covalent-length
parse/alt-conf artifacts) are never treated as H-bonds.

Shared by:
  - scripts/reference/build_prosco_distributions.py  (distance ProSco table)
  - scripts/reference/build_z_tables.py               (distance Z' table)
  - scripts/reference/calculate_penalty_weights.py    (distance severity)
  - src/pair_finder/scoring/issues.py                 (runtime distance severity)

Keep build and runtime in sync by routing both through this module.
"""

from __future__ import annotations

from collections.abc import Iterable

from .canonical import canonicalize_bp

# Cells judged on ALL base-base bonds rather than the primary bond only.
WC_DISTANCE_CELLS: frozenset[tuple[str, str]] = frozenset({
    ("C-G", "cWW"), ("A-U", "cWW"),
})

# "Strong" base-base bond cutoff — mirrors max_hbond_distance in the finder
# config (pair_finder.finder.finder.FinderConfig).
MAX_STRONG_DISTANCE = 3.5

# Below this a heavy-atom contact is not a real H-bond (covalent-length
# parse/alt-conf artifact); excluded everywhere.
MIN_PHYSICAL_DISTANCE = 2.2

# Base-ring atoms — mirrors pair_finder.scoring.features._BASE_ATOMS so the
# distance distribution is built from the same bonds the scorer evaluates.
BASE_ATOMS: frozenset[str] = frozenset({
    "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6",
})


def is_wc_distance_cell(bp_type: str | None, lw_class: str | None) -> bool:
    """True if (canonical bp_type, lw_class) is a Watson-Crick distance cell."""
    cb, _ = canonicalize_bp(bp_type, lw_class)
    return (cb, lw_class) in WC_DISTANCE_CELLS


def distances_to_score(distances: Iterable[float | None],
                       bp_type: str | None, lw_class: str | None) -> list[float]:
    """Runtime: the base-base distance(s) whose severity should be taken (max).

    WC cells return every physical bond; other cells return the single
    strongest bond within [MIN_PHYSICAL_DISTANCE, MAX_STRONG_DISTANCE], or an
    empty list when the pair has no strong bond.
    """
    phys = [d for d in distances if d is not None and d >= MIN_PHYSICAL_DISTANCE]
    if not phys:
        return []
    if is_wc_distance_cell(bp_type, lw_class):
        return phys
    strong = [d for d in phys if d <= MAX_STRONG_DISTANCE]
    return [min(strong)] if strong else []


def select_distance_rows(hb_df):
    """Build side: the subset of hbond rows that define the distance distribution.

    Requires columns: pdb_id, pair_idx, bp_type, lw_class, donor_atom,
    acceptor_atom, distance. Base-atom filtered with a physical floor, then
    WC cells keep all bonds while other cells keep only the strongest bond
    (<= MAX_STRONG_DISTANCE) per pair. bp_type is returned unchanged (callers
    canonicalize at grouping time).
    """
    import pandas as pd

    df = hb_df[
        hb_df["donor_atom"].isin(BASE_ATOMS)
        & hb_df["acceptor_atom"].isin(BASE_ATOMS)
        & (hb_df["distance"] >= MIN_PHYSICAL_DISTANCE)
    ].copy()
    canon_bp = [canonicalize_bp(b, lw)[0] for b, lw in zip(df["bp_type"], df["lw_class"])]
    is_wc = pd.Series(
        [(cb, lw) in WC_DISTANCE_CELLS for cb, lw in zip(canon_bp, df["lw_class"])],
        index=df.index,
    )
    wc_rows = df[is_wc]
    other = df[~is_wc]
    strong_idx = other.groupby(["pdb_id", "pair_idx"])["distance"].idxmin()
    strong = other.loc[strong_idx]
    strong = strong[strong["distance"] <= MAX_STRONG_DISTANCE]
    return pd.concat([wc_rows, strong])


def base_base_rows(hb_df):
    """Build side: base-base bonds only (for the angle distributions, to match
    the scorer which evaluates angles on base-base bonds)."""
    return hb_df[
        hb_df["donor_atom"].isin(BASE_ATOMS)
        & hb_df["acceptor_atom"].isin(BASE_ATOMS)
    ]
