"""Supplemental paired selection audit with comparator-conditioned candidate credit.

Each comparison retains FastRNA's selected output and adds classified discarded
candidates only when that comparator reports the same residue pair. FastRNA's own
LW label is preserved, including disagreements. Eligible candidate dispositions
are selected, competition, ineligible, stacking, and low_score. Geometrically
invalid candidates and adjacent residues are outside this benchmark.

This audit excludes mismatched candidate/output archives. Figure 2A retains its
archived classification metrics; Figure 2B-C and S3A-C retain the original full
pair-list cohort with selected-output fallback for those mismatched archives.
"""
from collections import Counter
from functools import cached_property, lru_cache
import csv
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

TOOLS = ["parse", "dssr", "fr3d", "rnaview", "barnaba", "rnapolis", "bpnet"]
LABELS = ["FastRNA", "DSSR", "FR3D", "RNAView", "barnaba", "RNApolis", "BPNET"]
IDENTIFIED = frozenset({"selected", "competition", "ineligible", "stacking", "low_score"})
NONCWW = ["tWW", "cWH", "tWH", "cWS", "tWS", "cHH", "tHH", "cHS", "tHS", "cSS", "tSS"]
CLASSES = ["cWW (canon)", "cWW (non-canon)", *NONCWW]
UNCLASSIFIED = "?"
_NUM = re.compile(r"^(-?\d+)([A-Za-z]?)$")


def residue(chain, number, icode=""):
    m = _NUM.fullmatch(str(number))
    if not chain or m is None or (m[2] and icode):
        raise ValueError(f"Uncanonicalizable residue ID: {chain!r} {number!r} {icode!r}")
    return str(chain), int(m[1]), icode or m[2]


def flip_label(label):
    def flip(lw):
        if lw.startswith("n"):
            return "n" + flip(lw[1:])
        return lw[0] + lw[2] + lw[1] if len(lw) == 3 else lw
    return "|".join(flip(x) for x in label.split("|"))


def pair_entry(a, b, label):
    return ((b, a), flip_label(label)) if b < a else ((a, b), label)


def adjacent(a, b):
    # Preserve the published benchmark's explicit sequence-neighbour exclusion.
    return a[0] == b[0] and abs(a[1] - b[1]) == 1


def is_lw(lw):
    return bool(lw) and len(lw) == 3 and lw[0] in "ct" and lw[1] in "WHS" and lw[2] in "WHS"


def classified(lw):
    return bool(lw) and all(is_lw(part) for part in lw.split("|"))


def single_family(labels):
    # Partial, ambiguous, and multiple labels are not definite class calls.
    return next(iter(labels)) if len(labels) == 1 and is_lw(next(iter(labels))) else None


def canon_lw(lw):
    order = {"W": 0, "H": 1, "S": 2}
    return lw[0] + lw[2] + lw[1] if order[lw[1]] > order[lw[2]] else lw


def ingest(rows):
    pairs = {}
    for a, b, lw in rows:
        if adjacent(a, b):
            continue
        key, label = pair_entry(a, b, lw or UNCLASSIFIED)
        pairs.setdefault(key, set()).add(label)
    return pairs


def single_family_units(pairsets):
    return {(pdb, key): lw for pdb, pairs in pairsets.items()
            for key, labels in pairs.items() if (lw := single_family(labels)) is not None}


def safe_ratio(n, d):
    return n / d if d else float("nan")


class Comparison:
    def __init__(self, data_dir):
        self.data_dir = Path(data_dir)

    def read(self, name):
        return json.loads((self.data_dir / "comparison" / name).read_text())

    @cached_property
    def candidates(self):
        raw = self.read("parse_candidates.json")
        if raw.get("_meta", {}).get("label_order") != "input_residue_order":
            raise ValueError("Regenerate parse_candidates.json with the orientation-corrected "
                             "dump-selection export; see scripts/validation/extract_parse_candidates.py")
        return raw

    @cached_property
    def excluded(self):
        return frozenset(self.candidates["_meta"].get("fallback_to_output", []))

    @cached_property
    def cohort(self):
        selected = self.read("parse_pairs.json")
        return frozenset(pdb for pdb, record in selected.items()
                         if record.get("status") == "ok" and pdb not in self.excluded
                         and self.candidates.get(pdb, {}).get("status") == "ok")

    @lru_cache(maxsize=None)
    def tool(self, name):
        if name == "fr3d":
            out = {}
            for path in sorted((self.data_dir / "comparison/fr3d").glob("[0-9A-Z]*.csv")):
                if path.stem not in self.cohort:
                    continue
                rows = []
                with path.open() as fh:
                    for row in csv.reader(fh):
                        if len(row) != 3:
                            continue
                        a, b = row[0].split("|"), row[2].split("|")
                        if min(len(a), len(b)) < 5 or a[1] != "1" or b[1] != "1":
                            continue
                        rows.append((residue(a[2], a[4], a[7] if len(a) > 7 else ""),
                                     residue(b[2], b[4], b[7] if len(b) > 7 else ""), row[1]))
                out[path.stem] = ingest(rows)
            return out
        raw = self.read(f"{name}_pairs.json")
        return {pdb: ingest((residue(r[0], r[1]), residue(r[2], r[3]), r[4])
                            for r in rec["pairs"])
                for pdb, rec in raw.items() if pdb in self.cohort and rec.get("status") == "ok"}

    @cached_property
    def selected(self):
        return {pdb: {k: labels for k, labels in pairs.items()
                      if all(classified(lw) for lw in labels)}
                for pdb, pairs in self.tool("parse").items()}

    @cached_property
    def identified(self):
        out = {}
        for pdb in sorted(self.cohort):
            rows = self.candidates[pdb]["pairs"]
            pairs = ingest((residue(a, n), residue(b, m), lw)
                           for a, n, b, m, lw, disposition in rows
                           if disposition in IDENTIFIED and classified(lw))
            # The archived final output is authoritative on already-selected pairs.
            # This also preserves its explicit ambiguity labels. No comparator is used.
            pairs.update(self.selected[pdb])
            out[pdb] = pairs
        return out

    def pairsets(self):
        return {t: self.identified if t == "parse" else self.tool(t) for t in TOOLS}

    @lru_cache(maxsize=None)
    def credited(self, comparator):
        """Selected output plus discarded candidates reported by this comparator."""
        reference = self.tool(comparator)
        return {pdb: {**{key: labels for key, labels in self.identified[pdb].items()
                         if key in reference.get(pdb, {})}, **selected}
                for pdb, selected in self.selected.items()}

    def pair_lists(self, a, b, stage="credited"):
        if stage not in {"credited", "selected", "identified"}:
            raise ValueError(f"Unknown comparison stage: {stage}")

        def pairs(tool, partner):
            if tool != "parse":
                return self.tool(tool)
            return self.credited(partner) if stage == "credited" else getattr(self, stage)

        return pairs(a, b), pairs(b, a)

    @cached_property
    def residue_names(self):
        return self.read("residue_names.json")

    def canonical(self, pdb, key, labels):
        names = self.residue_names.get(pdb, {})
        bases = [names.get(f"{c}|{n}|{i}") for c, n, i in key]
        return (None not in bases and frozenset(bases) in {frozenset("GC"), frozenset("AU")}
                and any("cWW" in lw.split("|") for lw in labels))

    def class_group(self, pdb, key, lw):
        if lw == "cWW":
            return "cWW (canon)" if self.canonical(pdb, key, {lw}) else "cWW (non-canon)"
        return canon_lw(lw)

    @lru_cache(maxsize=None)
    def detection_counts(self, a="parse", b="dssr", stage="credited"):
        """Per-structure intersection/union counts, canonical then noncanonical."""
        pa, pb = self.pair_lists(a, b, stage)
        pdbs = sorted(pa.keys() & pb.keys())
        shared, union = np.zeros((len(pdbs), 2)), np.zeros((len(pdbs), 2))
        excluded = 0
        for i, pdb in enumerate(pdbs):
            x, y = pa[pdb], pb[pdb]
            for key in x.keys() | y.keys():
                labels = x.get(key, set()) | y.get(key, set())
                if any("|" in lw for lw in labels):
                    excluded += 1
                    continue
                k = 0 if self.canonical(pdb, key, labels) else 1
                union[i, k] += 1
                shared[i, k] += key in x and key in y
        return pdbs, shared, union, excluded

    @lru_cache(maxsize=None)
    def class_counts(self, a="parse", b="dssr", stage="credited"):
        pa, pb = self.pair_lists(a, b, stage)
        pdbs = sorted(pa.keys() & pb.keys())
        same, shared = np.zeros((len(pdbs), 2)), np.zeros((len(pdbs), 2))
        for i, pdb in enumerate(pdbs):
            for key in pa[pdb].keys() & pb[pdb].keys():
                x, y = single_family(pa[pdb][key]), single_family(pb[pdb][key])
                if x is None or y is None:
                    continue
                k = 0 if self.canonical(pdb, key, {x, y}) else 1
                shared[i, k] += 1
                same[i, k] += x == y
        return pdbs, same, shared

    @cached_property
    def dssr_stage_counts(self):
        """Class-aware TP/predicted/reference counts, with one mask for all stages.

        On the pre-selection/DSSR union, a pair with an unresolved class in either
        present annotation is excluded from all stages. Absent pairs are not
        exclusions. Incorrect definite classes count against precision and recall.
        The reference denominator is fixed for identified, credited, and selected sets.
        """
        dssr = self.tool("dssr")
        rows, excluded = [], 0
        for pdb in sorted(self.identified.keys() & dssr.keys()):
            before, reference = self.identified[pdb], dssr[pdb]
            mask = {k for k in before.keys() | reference.keys()
                    if (k in before and single_family(before[k]) is None)
                    or (k in reference and single_family(reference[k]) is None)}
            excluded += len(mask)
            for stage, predicted in (("identified", before),
                                     ("credited", self.credited("dssr")[pdb]),
                                     ("selected", self.selected[pdb])):
                tp, pred, ref = Counter(), Counter(), Counter()
                for key in predicted.keys() | reference.keys():
                    if key in mask:
                        continue
                    x = single_family(predicted[key]) if key in predicted else None
                    y = single_family(reference[key]) if key in reference else None
                    if x is not None:
                        pred[self.class_group(pdb, key, x)] += 1
                    if y is not None:
                        ref[self.class_group(pdb, key, y)] += 1
                    if x is not None and x == y:
                        tp[self.class_group(pdb, key, x)] += 1
                for family in CLASSES:
                    rows.append((pdb, stage, family, tp[family], pred[family], ref[family]))
        return pd.DataFrame(rows, columns=["pdb", "stage", "family", "tp", "predicted", "reference"]), excluded

    def metrics(self, stage="credited"):
        df, _ = self.dssr_stage_counts
        totals = df[df.stage == stage].groupby("family")[["tp", "predicted", "reference"]].sum()
        return {family: {"recall_cls": safe_ratio(r.tp, r.reference),
                         "prec_cls": safe_ratio(r.tp, r.predicted),
                         "total_dssr": int(r.reference), "total_python": int(r.predicted),
                         "matched": int(r.tp)} for family, r in totals.iterrows()}

    def selection_effect(self):
        df, _ = self.dssr_stage_counts
        df = df.assign(kind=np.where(df.family == "cWW (canon)", "canonical", "non-canonical"))
        out = df.groupby(["stage", "kind"])[["tp", "predicted", "reference"]].sum()
        out["recall"] = out.tp / out.reference.replace(0, np.nan)
        out["precision"] = out.tp / out.predicted.replace(0, np.nan)
        return out

    def disposition_counts(self):
        counts = Counter()
        for pdb in self.cohort:
            seen = set()
            for a, n, b, m, lw, disposition in self.candidates[pdb]["pairs"]:
                if disposition not in IDENTIFIED or not classified(lw):
                    continue
                r1, r2 = residue(a, n), residue(b, m)
                if adjacent(r1, r2):
                    continue
                key, _ = pair_entry(r1, r2, lw)
                if key in seen:
                    raise ValueError(f"Duplicate candidate {pdb} {key}")
                seen.add(key)
                counts["selected" if key in self.selected[pdb] else disposition] += 1
        return pd.Series({k: counts[k] for k in sorted(IDENTIFIED)}, name="classified_pairs")
