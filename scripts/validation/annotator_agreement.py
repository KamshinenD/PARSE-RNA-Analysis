"""Inter-annotator agreement across PARSE and the five-referee panel.

Raters (6): PARSE (data/comparison/parse_pairs.json, run_parse_annotator.py)
plus the five voting referees — published FR3D (BGSU Hub CSVs), RNAView,
barnaba, RNApolis, BPNET. fr3d_local is non-voting and excluded. DSSR is NOT
a rater here: its full per-structure pair lists are not on disk, and its
appearances inside the data/comparison/{class}.json `differences` records
are a disagreement-only (biased) sample, so it is omitted from every matrix
below and this is stated in the output JSON.

All loading goes through parse_analysis.annotators (the canonical-id layer):
units are order-normalized PairKeys of ResidueId (multi-char chains and
insertion codes intact), labels expressed in the canonical key direction.

Metrics, written to data/validation/annotator_agreement.json:
  1. detection_jaccard: per rater pair, pooled (micro) Jaccard of detected
     pair-key sets over structures both cover with status ok, plus the
     median (proper even-n median) per-structure Jaccard.
  2. classification: per rater pair, over pair-keys both report with a single
     LW family (null/partial labels skipped): percent agreement + Cohen's
     kappa. Two variants: "strict" excludes PARSE-ambiguous pairs;
     "lenient" scores a PARSE-ambiguous "X|Y" pair as agreeing when the other
     rater's family is in {X, Y} (for kappa, PARSE's nominal label is then
     the matching candidate, else the first candidate).
  3. krippendorff_alpha (nominal) across all 6 raters: units = (structure,
     pair-key) labeled with a single family by >= 2 raters; PARSE-ambiguous
     pairs count as missing for PARSE. Leave-one-out alphas are PAIRED: for
     each dropped rater, the SAME unit set (units still >= 2 raters without
     the dropped one) is scored with and without that rater, and the paired
     delta reported. Alpha is CONDITIONAL on co-detection (units detected and
     single-family-labeled by >= 2 raters); the JSON carries that fraction
     and a note, plus per-rater label-loss (pairs dropped for lacking a
     unique single-family label).
  4. parse_consensus_support: for each PARSE-detected pair, k = number of
     the 5 referees covering that structure that also detect it, split by
     PARSE-unambiguous vs PARSE-ambiguous.

Run AFTER run_parse_annotator.py (and the referee runners):
    python scripts/validation/annotator_agreement.py
"""

import doctest
import json
import sys
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from statistics import median

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))  # vendored annotators
from _vendor import annotators  # noqa: E402

OUT = REPO / "data" / "validation" / "annotator_agreement.json"
RATERS = ("parse", "fr3d", "rnaview", "barnaba", "rnapolis", "bpnet")
REFEREES = RATERS[1:]
DSSR_NOTE = (
    "DSSR is not a rater: full DSSR pair lists are not on disk; DSSR appears only inside "
    "the comparison `differences` records, a disagreement-only (biased) sample, so it is "
    "excluded from all matrices. The primary panel is the 6 raters listed."
)
ALPHA_NOTE = (
    "Alpha and the kappas are CONDITIONAL on co-detection: they score only (structure, "
    "pair-key) units that >= 2 raters both detected and gave a unique single-family label "
    "(co_detection_fraction of the units any rater labeled). Detection disagreement is "
    "reported separately in detection_jaccard and is NOT inside alpha."
)
COVERAGE_NOTE = (
    "Coverage is whatever each tool's *_pairs.json reports as status ok; RNAView reaches "
    "many-chain structures through single-char chain aliases and, where needed, chain "
    "groups, and barnaba reaches mmCIF-only structures through synthetic PDB files "
    "(run_third_annotators.py). RNAView and barnaba are icode-blind (their output formats "
    "drop insertion codes), so pairs on insertion-coded residues cannot match them and "
    "count against their detection rates."
)


def load_raters() -> dict[str, dict[str, annotators.PairSets]]:
    """All 6 raters as {rater: {pdb: canonical PairSets}} (ok structures only)."""
    raters: dict[str, dict[str, annotators.PairSets]] = {"parse": annotators.parse_pairs()}
    raters["fr3d"] = annotators.fr3d_pairs()
    for tool in REFEREES[1:]:
        raters[tool] = annotators.tool_pairs(tool)
    return raters


def single_family(labels: set[str]) -> str | None:
    """The unique single-family LW label, or None (ambiguous/partial/multi)."""
    if any("|" in lw for lw in labels):
        return None
    fams = {lw for lw in labels if len(lw) == 3 and lw[0] in "ct"}
    return next(iter(fams)) if len(fams) == 1 else None


def ambiguous_candidates(labels: set[str]) -> set[str] | None:
    """The candidate families of a PARSE-ambiguous label set, or None."""
    for lw in labels:
        if "|" in lw:
            return set(lw.split("|"))
    return None


def _fmt(x) -> str:
    return "  n/a" if x is None else f"{x:.3f}"


def detection_matrix(raters: dict) -> dict:
    """Pooled + median per-structure Jaccard of detected pair-key sets."""
    jac = {a: {a: 1.0} for a in RATERS}
    med = {a: {a: 1.0} for a in RATERS}
    n_shared = {a: {a: len(raters[a])} for a in RATERS}
    for a, b in combinations(RATERS, 2):
        shared = raters[a].keys() & raters[b].keys()
        inter = union = 0
        per_structure = []
        for p in shared:
            i = len(raters[a][p].keys() & raters[b][p].keys())
            u = len(raters[a][p].keys() | raters[b][p].keys())
            inter, union = inter + i, union + u
            if u:
                per_structure.append(i / u)
        jac[a][b] = jac[b][a] = round(inter / union, 4) if union else None
        mid = median(per_structure) if per_structure else None
        med[a][b] = med[b][a] = round(mid, 4) if mid is not None else None
        n_shared[a][b] = n_shared[b][a] = len(shared)
    return {
        "jaccard_pooled": jac,
        "jaccard_median_per_structure": med,
        "n_shared_structures": n_shared,
        "note_coverage": COVERAGE_NOTE,
    }


def cohen_kappa(pairs: list[tuple[str, str]]) -> float | None:
    """Cohen's kappa for two raters' nominal labels (paired list)."""
    n = len(pairs)
    if n == 0:
        return None
    po = sum(1 for x, y in pairs if x == y) / n
    m1 = Counter(x for x, _ in pairs)
    m2 = Counter(y for _, y in pairs)
    pe = sum(m1[c] * m2.get(c, 0) for c in m1) / (n * n)
    return 1.0 if pe == 1.0 else (po - pe) / (1 - pe)


def _labelled_units(ra: dict, rb: dict, a_is_parse: bool, lenient: bool) -> list:
    """Label pairs [(fa, fb), ...] for pair-keys both raters report."""
    out = []
    for key in ra.keys() & rb.keys():
        fb = single_family(rb[key])
        if fb is None:
            continue
        fa = single_family(ra[key])
        if fa is None and lenient and a_is_parse:
            cands = ambiguous_candidates(ra[key])
            if cands:
                fa = fb if fb in cands else min(cands)
        if fa is not None:
            out.append((fa, fb))
    return out


def classification_pair(raters: dict, a: str, b: str, lenient: bool) -> dict:
    """Percent agreement + kappa for one rater pair, one variant."""
    pairs: list[tuple[str, str]] = []
    for pdb in raters[a].keys() & raters[b].keys():
        pairs.extend(_labelled_units(raters[a][pdb], raters[b][pdb], a == "parse", lenient))
    n = len(pairs)
    agree = sum(1 for x, y in pairs if x == y)
    kappa = cohen_kappa(pairs)
    return {
        "n_pairs": n,
        "percent_agreement": round(100 * agree / n, 2) if n else None,
        "kappa": round(kappa, 4) if kappa is not None else None,
    }


def classification_matrices(raters: dict) -> dict:
    """Both classification variants for every rater pair."""
    out: dict = {"strict": {}, "lenient": {}}
    for variant, lenient in (("strict", False), ("lenient", True)):
        cells = {a: {} for a in RATERS}
        for a, b in combinations(RATERS, 2):
            cell = classification_pair(raters, a, b, lenient)
            cells[a][b] = cells[b][a] = cell
        out[variant] = cells
    return out


def krippendorff_alpha(units: list[list[str]]) -> float:
    """Krippendorff's alpha, nominal metric, from per-unit value lists.

    Units with fewer than 2 values are ignored (missing data allowed).

    >>> round(krippendorff_alpha([["a", "a"], ["a", "b"], ["b", "b"]]), 6)
    0.444444
    >>> krippendorff_alpha([["a", "a"], ["b", "b", "b"]])
    1.0
    """
    coincidence: Counter = Counter()
    for values in units:
        m = len(values)
        if m < 2:
            continue
        for i, v in enumerate(values):
            for j, w in enumerate(values):
                if i != j:
                    coincidence[(v, w)] += 1.0 / (m - 1)
    totals: Counter = Counter()
    for (v, _w), x in coincidence.items():
        totals[v] += x
    n = sum(totals.values())
    d_o = sum(x for (v, w), x in coincidence.items() if v != w) / n
    d_e = sum(totals[v] * totals[w] for v in totals for w in totals if v != w) / (n * (n - 1))
    return 1.0 if d_e == 0 else 1 - d_o / d_e


def alpha_units(raters: dict) -> dict[tuple, dict[str, str]]:
    """(pdb, pair-key) -> {rater: family} for single-family labels only."""
    units: dict[tuple, dict[str, str]] = defaultdict(dict)
    for rater in RATERS:
        for pdb, pairs in raters[rater].items():
            for key, labels in pairs.items():
                fam = single_family(labels)
                if fam is not None:
                    units[(pdb, key)][rater] = fam
    return units


def label_loss(raters: dict) -> dict:
    """Per-rater count of detected pairs dropped for lacking a single family."""
    out = {}
    for rater in RATERS:
        n_pairs = n_labelled = 0
        for pairs in raters[rater].values():
            n_pairs += len(pairs)
            n_labelled += sum(1 for labels in pairs.values() if single_family(labels))
        out[rater] = {
            "n_detected_pairs": n_pairs,
            "n_single_family": n_labelled,
            "n_dropped": n_pairs - n_labelled,
            "frac_dropped": round((n_pairs - n_labelled) / n_pairs, 4) if n_pairs else None,
        }
    return out


def alpha_report(raters: dict) -> dict:
    """Alpha for all 6 raters plus PAIRED leave-one-out alphas."""
    units = alpha_units(raters)
    all_values = [list(v.values()) for v in units.values() if len(v) >= 2]
    report = {
        "note": ALPHA_NOTE,
        "n_units_ge1_rater": len(units),
        "n_units_ge2_raters": len(all_values),
        "co_detection_fraction": round(len(all_values) / len(units), 4) if units else None,
        "alpha_all_6": round(krippendorff_alpha(all_values), 4),
        "leave_one_out_paired": {},
    }
    for drop in RATERS:
        paired = [v for v in units.values()
                  if sum(1 for r in v if r != drop) >= 2]
        alpha_with = krippendorff_alpha([list(v.values()) for v in paired])
        alpha_without = krippendorff_alpha(
            [[f for r, f in v.items() if r != drop] for v in paired])
        report["leave_one_out_paired"][f"drop_{drop}"] = {
            "n_units_paired": len(paired),
            "alpha_with": round(alpha_with, 4),
            "alpha_without": round(alpha_without, 4),
            "paired_delta_without_minus_with": round(alpha_without - alpha_with, 4),
        }
    return report


def consensus_support(raters: dict) -> dict:
    """k-of-panel referee detection distribution for PARSE pairs."""
    dists = {"unambiguous": Counter(), "ambiguous": Counter()}
    panel5 = {"unambiguous": Counter(), "ambiguous": Counter()}
    for pdb, pairs in raters["parse"].items():
        panel = [r for r in REFEREES if pdb in raters[r]]
        for key, labels in pairs.items():
            k = sum(1 for r in panel if key in raters[r][pdb])
            group = "ambiguous" if ambiguous_candidates(labels) else "unambiguous"
            dists[group][(len(panel), k)] += 1
            if len(panel) == 5:
                panel5[group][k] += 1
    return {
        "k_by_panel": {g: {f"panel_{p}_k_{k}": n for (p, k), n in sorted(c.items())}
                       for g, c in dists.items()},
        "k_panel5_only": {g: {f"k_{k}": n for k, n in sorted(c.items())}
                          for g, c in panel5.items()},
    }


def sanity_checks() -> dict:
    """Dedup record count + the consensus majority the normalization must match."""
    n_records = sum(1 for _ in annotators.disagreement_records())
    with open(REPO / "data" / "validation" / "annotator_consensus.json") as fh:
        consensus = json.load(fh)
    return {
        "deduped_disagreement_records": n_records,
        "expected_deduped": 123022,
        "consensus_panel5_majority": consensus["class_disagreement_majority"]["panel_5"],
    }


def print_summary(result: dict) -> None:
    for kind in ("jaccard_pooled", "jaccard_median_per_structure"):
        det = result["detection_jaccard"][kind]
        print(f"\ndetection {kind} (shared-ok structures):")
        for a in RATERS:
            print("  " + a.ljust(9) + " ".join(
                f"{b}={_fmt(det[a][b])}" for b in RATERS if b != a))
    for variant in ("strict", "lenient"):
        cells = result["classification"][variant]
        print(f"\nclassification, {variant} (pct agree / kappa / n):")
        for b in REFEREES:
            c = cells["parse"][b]
            print(f"  parse-{b:<9} {c['percent_agreement']}% / {c['kappa']} / n={c['n_pairs']}")
    alpha = result["krippendorff_alpha"]
    print(f"\nKrippendorff alpha (nominal, all 6): {alpha['alpha_all_6']} "
          f"over {alpha['n_units_ge2_raters']} units "
          f"(co-detection fraction {alpha['co_detection_fraction']})")
    for drop, cell in alpha["leave_one_out_paired"].items():
        print(f"  {drop:<14} with={cell['alpha_with']} without={cell['alpha_without']} "
              f"paired_delta={cell['paired_delta_without_minus_with']:+.4f} "
              f"(n={cell['n_units_paired']})")
    print("\nper-rater label loss (detected pairs without a unique family):")
    for rater, cell in result["label_loss"].items():
        print(f"  {rater:<9} dropped {cell['n_dropped']:>7} of {cell['n_detected_pairs']:>8} "
              f"({cell['frac_dropped']})")
    print("\nPARSE pair referee support (panel of 5 only):")
    for group, counts in result["parse_consensus_support"]["k_panel5_only"].items():
        tot = sum(counts.values())
        line = ", ".join(f"{k}: {100 * n / tot:.1f}%" for k, n in counts.items())
        print(f"  {group} (n={tot}): {line}")
    print(f"\nsanity: {result['sanity']}")


def main() -> None:
    failures = doctest.testmod(sys.modules[__name__]).failed
    if failures:
        raise SystemExit(f"{failures} doctest failure(s)")
    raters = load_raters()
    print({r: len(v) for r, v in raters.items()}, "structures with ok status")
    result = {
        "note_dssr": DSSR_NOTE,
        "raters": list(RATERS),
        "detection_jaccard": detection_matrix(raters),
        "classification": classification_matrices(raters),
        "krippendorff_alpha": alpha_report(raters),
        "label_loss": label_loss(raters),
        "parse_consensus_support": consensus_support(raters),
        "sanity": sanity_checks(),
    }
    OUT.write_text(json.dumps(result, indent=1))
    print_summary(result)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
