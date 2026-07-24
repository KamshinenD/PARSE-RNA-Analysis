"""Richardson 2008 RNA backbone suite classification.

Assigns a 7D backbone suite to one of 54 Richardson conformers (e.g. ``1a``,
``1g``, ``&a``) and computes a 0-1 ``suiteness`` score that measures how
close the suite sits to its conformer cluster center. Outliers that fail
triage / sieve / 4D / 7D screening receive ``conformer = "!!"`` and
``suiteness = 0``.

The 7D suite for residue *i* is::

    [delta(i-1), epsilon(i-1), zeta(i-1), alpha(i), beta(i), gamma(i), delta(i)]

All angles are expected in the 0-360 range (apply ``angle % 360`` before
calling).

Algorithm (from Richardson et al. 2008):
    1. Triage: epsilon, alpha, beta, zeta must lie in plausible ranges.
    2. Sieve: classify both deltas (C2'-endo / C3'-endo) and gamma (p/t/m).
       The (delta_prev, delta, gamma) triple selects a conformer bin.
    3. 4D screen: find the conformer in the bin whose [eps, zeta, alpha,
       beta] center is closest by L3 (Minkowski p=3) hyperellipsoid distance
       using width-normalized differences. Reject if d4 >= 1.
       Widths are conformer-specific:
         * satellites use sat_widths over the satellite baseline,
         * dominants in a bin containing satellites use the dom_widths from
           the satellite that gives the dominant the best fit,
         * ordinary conformers use the normal baseline.
    4. 7D test: same L3 distance, with the same conformer-specific widths,
       across all 7 angles. Reject if d7 >= 1.
    5. Suiteness: raised cosine of d7  →  ``(cos(pi*d7) + 1) / 2``,
       floored at 0.01.

Ported from score_structure_motifs/scoring/scorer.py and refactored to a
standalone, picklable classifier.

Reference: Richardson JS et al. (2008) RNA 14:465.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

DEFAULT_SUITES_PATH = Path(__file__).resolve().parent / "richardson_suites.json"

# Outlier marker (matches Richardson convention).
OUTLIER_CONFORMER = "!!"

# Names of the seven suite dihedrals, in the order the suite array carries them:
# [delta(i-1), epsilon(i-1), zeta(i-1), alpha(i), beta(i), gamma(i), delta(i)].
SUITE_ANGLE_NAMES = (
    "delta_prev", "epsilon_prev", "zeta_prev",
    "alpha", "beta", "gamma", "delta",
)


@dataclass
class SuiteResult:
    """Result of classifying one 7D backbone suite.

    Attributes:
        conformer: Name of the matched Richardson conformer (e.g. ``"1a"``)
            or ``"!!"`` for outliers.
        suiteness: 0-1 score, raised cosine of 7D L3 distance to cluster
            center. Floored at 0.01 for in-cluster points; 0 for outliers.
        bin: Conformer bin key formed from (delta_prev, delta, gamma)
            classifications, e.g. ``"33p"``. None when triage / sieve fail.
        is_outlier: True when the suite did not pass triage, sieve, 4D,
            or 7D screening.
        distance: 7D L3 distance to the matched cluster center, or None
            if no match was made.
        issue: Human-readable reason for outlier classification, or None
            when a match succeeded.
    """

    conformer: str
    suiteness: float
    bin: str | None
    is_outlier: bool
    distance: float | None
    issue: str | None


def _angle_diff(a: float, b: float) -> float:
    """Angular difference in 0-360 space, taking wraparound into account."""
    d = abs(a - b)
    return min(d, 360.0 - d)


def signed_angle_gap(value: float, center: float) -> float:
    """Signed minimal angular gap ``value - center`` in (-180, 180].

    Positive means ``value`` sits above ``center`` on the shorter arc, so a
    fix would *decrease* the angle; negative means the reverse. Used to tell a
    refiner which direction to rotate a deviating suite torsion.
    """
    return (value - center + 180.0) % 360.0 - 180.0


def _l3_distance(suite: np.ndarray, center: np.ndarray, widths: np.ndarray,
                 indices: range) -> float:
    """L3 (Minkowski p=3) hyperellipsoid distance over the given indices."""
    total = 0.0
    for k in indices:
        d = _angle_diff(suite[k], center[k])
        total += (d / widths[k]) ** 3
    return total ** (1.0 / 3.0)


class RichardsonClassifier:
    """Loads the Richardson conformer table once and reuses it.

    Construct one classifier per process and call ``classify(suite_angles)``
    for each residue. Loading the JSON file inside a tight loop is wasteful;
    the conformer index is also pre-built here so per-call lookups are O(1)
    in the bin and linear in the (small) number of conformers per bin.
    """

    def __init__(self, suites_path: Path | str = DEFAULT_SUITES_PATH):
        path = Path(suites_path)
        if not path.exists():
            raise FileNotFoundError(f"Richardson suites JSON not found: {path}")

        with open(path) as f:
            self._raw = json.load(f)

        # Pre-index conformers by their (deltaMinus, delta, gamma) bin
        # for fast lookup during sieve.
        self._by_bin: dict[str, list[dict]] = {}
        for conf in self._raw["conformers"]:
            self._by_bin.setdefault(conf["bin"], []).append(conf)

        self._normal_widths = np.array(self._raw["widths"]["normal"], dtype=float)
        self._satellite_widths = np.array(
            self._raw["widths"].get("satellite", self._raw["widths"]["normal"]),
            dtype=float,
        )
        self._sat_overrides: dict[str, dict] = self._raw.get("satellite_overrides", {})
        self._triage = self._raw["triage"]
        self._sieve = self._raw["sieve"]

        # Conformer name -> 7D cluster center, for per-angle deviation reporting.
        self._centers: dict[str, np.ndarray] = {
            conf["name"]: np.array(conf["angles"], dtype=float)
            for conf in self._raw["conformers"]
        }

    # ----------------------- deviation-reporting API ----------------------- #

    @property
    def normal_widths(self) -> np.ndarray:
        """Baseline per-axis tolerance widths (used to flag deviating angles)."""
        return self._normal_widths

    def center_for(self, conformer: str) -> np.ndarray | None:
        """7D cluster center for a conformer name, or None if unknown."""
        return self._centers.get(conformer)

    def nearest_conformer(self, suite_angles: np.ndarray | list[float]) -> str | None:
        """Name of the conformer whose center is closest to ``suite_angles``.

        Uses summed absolute angular gap (L1) over all seven axes. Used to give
        an outlier suite (no assigned conformer) a fix target: the rotamer it
        sits nearest to, even when it is too far to be classified into it.
        """
        suite = np.asarray(suite_angles, dtype=float)
        best_name, best_d = None, float("inf")
        for name, center in self._centers.items():
            d = sum(_angle_diff(suite[k], center[k]) for k in range(7))
            if d < best_d:
                best_d, best_name = d, name
        return best_name

    # ----------------------------- helpers ----------------------------- #

    def _classify_pucker(self, delta: float) -> int:
        """Classify sugar pucker: 3 = C3'-endo, 2 = C2'-endo, 0 = unclassified."""
        c3 = self._sieve["delta"]["C3_endo"]
        c2 = self._sieve["delta"]["C2_endo"]
        if c3[0] <= delta <= c3[1]:
            return 3
        if c2[0] <= delta <= c2[1]:
            return 2
        return 0

    def pucker(self, delta: float | None) -> int | None:
        """Public sugar-pucker classification from a residue's δ torsion.

        Returns 3 (C3'-endo), 2 (C2'-endo), 0 (between the two bands), or None
        when δ is unavailable. Uses the same δ sieve as suite classification
        (C3'-endo [60,105], C2'-endo [125,165]).
        """
        if delta is None:
            return None
        return self._classify_pucker(float(delta) % 360)

    def _classify_gamma(self, gamma: float) -> str:
        """Classify gamma rotamer: 'p' / 't' / 'm', or '' for unclassified."""
        for name, rng in self._sieve["gamma"].items():
            if rng[0] <= gamma <= rng[1]:
                return name
        return ""

    # ------------------ width resolution (dom/sat overrides) ------------------ #

    @staticmethod
    def _compose_widths(override: list[float], baseline: np.ndarray) -> np.ndarray:
        """Return baseline with non-zero override entries substituted in.

        Override widths use 0 to mean "fall back to baseline on this axis".
        Non-zero entries replace the baseline on the matching axis.
        """
        result = baseline.copy()
        for k, w in enumerate(override):
            if w > 0:
                result[k] = w
        return result

    def _widths_for_satellite(self, name: str) -> np.ndarray:
        """Effective widths for a satellite conformer (own sat_widths)."""
        ov = self._sat_overrides.get(name)
        if not ov:
            return self._satellite_widths
        return self._compose_widths(ov["sat_widths"], self._satellite_widths)

    def _widths_for_dominant(self, suite: np.ndarray,
                             dom_center: np.ndarray,
                             bin_confs: list[dict]) -> np.ndarray:
        """Effective widths for the dominant conformer in a bin.

        Each satellite in the bin defines its own dom_widths override
        describing how the dominant's territory should be reshaped when
        competing with *that* satellite. We try every override and keep
        the one yielding the smallest d4 — the dominant's best fit against
        any satellite-paired width set.
        """
        best_d = float("inf")
        best_widths: np.ndarray | None = None
        for sib in bin_confs:
            if sib.get("dominance") != "sat":
                continue
            ov = self._sat_overrides.get(sib["name"])
            if not ov:
                continue
            widths = self._compose_widths(ov["dom_widths"], self._normal_widths)
            d = _l3_distance(suite, dom_center, widths, range(1, 5))
            if d < best_d:
                best_d = d
                best_widths = widths
        return best_widths if best_widths is not None else self._normal_widths

    def _resolve_widths(self, conf: dict, bin_confs: list[dict],
                        suite: np.ndarray) -> np.ndarray:
        """Pick the effective width array for one conformer in its bin."""
        dominance = conf.get("dominance", "ord")
        if dominance == "sat":
            return self._widths_for_satellite(conf["name"])
        if dominance == "dom":
            center = np.asarray(conf["angles"], dtype=float)
            return self._widths_for_dominant(suite, center, bin_confs)
        return self._normal_widths

    # ----------------------------- public API ----------------------------- #

    def classify(self, suite_angles: np.ndarray | list[float]) -> SuiteResult:
        """Classify a 7D suite into a Richardson conformer.

        Args:
            suite_angles: 7-element array
                ``[delta_prev, epsilon_prev, zeta_prev, alpha, beta, gamma, delta]``
                with each angle in [0, 360). Caller is responsible for the
                ``angle % 360`` normalization.

        Returns:
            SuiteResult with conformer name, suiteness, and diagnostic info.
            On any rejection, ``conformer == "!!"``, ``suiteness == 0``,
            ``is_outlier == True``, and ``issue`` describes which check
            failed.
        """
        suite = np.asarray(suite_angles, dtype=float)
        if suite.shape != (7,):
            raise ValueError(
                f"suite_angles must have shape (7,), got {suite.shape}"
            )

        # ---- Triage: epsilon, alpha, beta, zeta in plausible ranges ----
        for idx, name in [(1, "epsilon"), (3, "alpha"), (4, "beta"), (2, "zeta")]:
            lo, hi = self._triage[name]
            if not (lo <= suite[idx] <= hi):
                return SuiteResult(
                    conformer=OUTLIER_CONFORMER, suiteness=0.0, bin=None,
                    is_outlier=True, distance=None,
                    issue=f"{name}={suite[idx]:.1f} outside [{lo}, {hi}]",
                )

        # ---- Sieve: classify both deltas and gamma ----
        pucker_dm = self._classify_pucker(suite[0])  # delta_prev
        pucker_d = self._classify_pucker(suite[6])   # delta
        gamma_class = self._classify_gamma(suite[5])

        if pucker_dm == 0 or pucker_d == 0 or gamma_class == "":
            return SuiteResult(
                conformer=OUTLIER_CONFORMER, suiteness=0.0, bin=None,
                is_outlier=True, distance=None,
                issue=(f"sieve fail: pucker_dm={pucker_dm} pucker_d={pucker_d} "
                       f"gamma='{gamma_class}'"),
            )

        bin_key = f"{pucker_dm}{pucker_d}{gamma_class}"
        conformers = self._by_bin.get(bin_key, [])
        if not conformers:
            return SuiteResult(
                conformer=OUTLIER_CONFORMER, suiteness=0.0, bin=bin_key,
                is_outlier=True, distance=None,
                issue=f"no conformers in bin {bin_key}",
            )

        # ---- 4D screen on (epsilon, zeta, alpha, beta) ----
        # Each conformer is scored against its conformer-specific widths
        # (satellite vs dominant vs ordinary). The widths chosen for the
        # winner are reused for the 7D test, keeping the metric consistent.
        best_d4 = float("inf")
        best_conf = None
        best_widths: np.ndarray | None = None
        for conf in conformers:
            center = np.array(conf["angles"])
            widths = self._resolve_widths(conf, conformers, suite)
            d4 = _l3_distance(suite, center, widths, range(1, 5))
            if d4 < best_d4:
                best_d4 = d4
                best_conf = conf
                best_widths = widths

        if best_conf is None or best_d4 >= 1.0:
            return SuiteResult(
                conformer=OUTLIER_CONFORMER, suiteness=0.0, bin=bin_key,
                is_outlier=True, distance=None,
                issue=f"no 4D match (best d4={best_d4:.2f})",
            )

        # ---- 7D test ----
        center = np.array(best_conf["angles"])
        d7 = _l3_distance(suite, center, best_widths, range(7))
        if d7 > 1.0:
            return SuiteResult(
                conformer=OUTLIER_CONFORMER, suiteness=0.0, bin=bin_key,
                is_outlier=True, distance=round(d7, 3),
                issue=f"7D reject {best_conf['name']} (d7={d7:.2f})",
            )

        # ---- Suiteness via raised cosine, floored at 0.01 ----
        suiteness = (math.cos(math.pi * d7) + 1.0) / 2.0
        if suiteness < 0.01:
            suiteness = 0.01

        return SuiteResult(
            conformer=best_conf["name"],
            suiteness=round(suiteness, 3),
            bin=bin_key,
            is_outlier=False,
            distance=round(d7, 3),
            issue=None,
        )


# Module-level helper for cheap chi classification (no Richardson data needed).

def classify_chi(chi: float | None) -> str | None:
    """Classify the glycosidic chi torsion as 'anti', 'syn', or None.

    Range matches X3DNA-DSSR convention:
        anti: chi <= -90 or chi >= +90
        syn:  -90 < chi < +90

    Args:
        chi: Glycosidic angle in degrees (-180, +180 range), or None when
            the dihedral could not be computed.

    Returns:
        'anti' | 'syn' for valid input; None when chi is None.
    """
    if chi is None:
        return None
    if chi <= -90.0 or chi >= 90.0:
        return "anti"
    return "syn"
