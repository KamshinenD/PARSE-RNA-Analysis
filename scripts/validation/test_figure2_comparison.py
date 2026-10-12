"""Small, hand-counted regression cases for the publication benchmark.

Run: python -m unittest discover -s scripts/validation -p 'test_figure2_comparison.py'
"""
import json
from pathlib import Path
import tempfile
import unittest

from figure2_comparison import Comparison, flip_label, ingest, residue


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        root = Path(self.tmp.name) / "comparison"
        root.mkdir()

        def write(name, value):
            (root / name).write_text(json.dumps(value))

        selected = [["A", "1", "A", "4", "cWW"], ["A", "40", "A", "30", "cSH"]]
        candidates = [
            [*selected[0], "selected"], [*selected[1], "selected"],
            ["A", "10", "A", "20", "tWH", "competition"],
            ["A", "50", "A", "60", "tWW", "ineligible"],  # unsupported prediction
            ["A", "70", "A", "80", "cWW|tWW", "stacking"],
            ["A", "90", "A", "100", None, "ineligible"],
            ["A", "110", "A", "120", "cWW", "invalid"],
        ]
        write("parse_pairs.json", {p: {"status": "ok", "pairs": selected} for p in ["P1", "OLD"]})
        write("parse_candidates.json", {
            "P1": {"status": "ok", "pairs": candidates},
            "OLD": {"status": "ok", "pairs": candidates},
            "_meta": {"label_order": "input_residue_order", "fallback_to_output": ["OLD"]},
        })
        write("dssr_pairs.json", {"P1": {"status": "ok", "pairs": [
            ["A", "1", "A", "4", "cWW"], ["A", "10", "A", "20", "tWH"],
            ["A", "30", "A", "40", "cWH"],  # detected, but class disagrees
            ["A", "70", "A", "80", "tWW"],  # same unresolved-label mask at both stages
        ]}})
        write("barnaba_pairs.json", {"P1": {"status": "ok", "pairs": []}})
        write("residue_names.json", {"P1": {"A|1|": "G", "A|4|": "C"}})
        self.c = Comparison(root.parent)

    def test_identification_does_not_depend_on_comparator(self):
        original = {k: set(v) for k, v in self.c.identified["P1"].items()}
        self.c.detection_counts("parse", "dssr")
        self.c.detection_counts("parse", "barnaba")
        self.assertEqual(self.c.identified["P1"], original)
        self.assertEqual(len(original), 5)
        # Unsupported tWW call remains, and costs precision.
        self.assertEqual(self.c.metrics("identified")["tWW"]["total_python"], 1)
        self.assertEqual(self.c.metrics("identified")["tWW"]["prec_cls"], 0)

    def test_credit_requires_the_current_comparator_to_report_the_pair(self):
        credited = self.c.credited("dssr")["P1"]
        supported = (residue("A", "10"), residue("A", "20"))
        unsupported = (residue("A", "50"), residue("A", "60"))
        self.assertIn(supported, credited)
        self.assertNotIn(unsupported, credited)
        self.assertEqual(len(credited), 4)
        self.assertEqual(self.c.credited("barnaba")["P1"], self.c.selected["P1"])
        self.assertEqual(self.c.metrics()["tWW"]["total_python"], 0)
        self.assertEqual(self.c.metrics()["tWH"]["matched"], 1)

    def test_credit_keeps_our_class_even_when_the_comparator_disagrees(self):
        key = (residue("A", "10"), residue("A", "20"))
        self.c.tool("dssr")["P1"][key] = {"cWH"}
        self.assertEqual(self.c.credited("dssr")["P1"][key], {"tWH"})
        self.assertEqual(self.c.metrics()["tWH"]["total_python"], 1)
        self.assertEqual(self.c.metrics()["tWH"]["matched"], 0)

    def test_selected_pairs_remain_even_without_comparator_support(self):
        _, shared, union, _ = self.c.detection_counts("parse", "barnaba")
        self.assertEqual(shared.sum(), 0)
        self.assertEqual(union.sum(), 2)

    def test_conditioned_detection_is_symmetric_and_class_counts_are_preserved(self):
        forward = self.c.detection_counts("parse", "dssr")
        reverse = self.c.detection_counts("dssr", "parse")
        self.assertEqual(forward[0], reverse[0])
        self.assertTrue((forward[1] == reverse[1]).all())
        self.assertTrue((forward[2] == reverse[2]).all())
        self.assertEqual(forward[1].sum(), 3)
        self.assertEqual(forward[2].sum(), 3)
        self.assertEqual(self.c.detection_counts(stage="identified")[2].sum(), 4)
        credited = self.c.class_counts()
        identified = self.c.class_counts(stage="identified")
        self.assertTrue((credited[1] == identified[1]).all())
        self.assertTrue((credited[2] == identified[2]).all())

    def test_same_reference_denominator_across_selection(self):
        effect = self.c.selection_effect()
        for kind in ["canonical", "non-canonical"]:
            self.assertEqual(effect.loc[("identified", kind), "reference"],
                             effect.loc[("selected", kind), "reference"])
            self.assertEqual(effect.loc[("credited", kind), "reference"],
                             effect.loc[("selected", kind), "reference"])
        self.assertAlmostEqual(effect.loc[("identified", "non-canonical"), "precision"], 1 / 3)
        self.assertAlmostEqual(effect.loc[("identified", "non-canonical"), "recall"], 1 / 2)
        self.assertEqual(effect.loc[("selected", "non-canonical"), "recall"], 0)
        self.assertAlmostEqual(effect.loc[("credited", "non-canonical"), "precision"], 1 / 2)
        self.assertAlmostEqual(effect.loc[("credited", "non-canonical"), "recall"], 1 / 2)
        self.assertEqual(self.c.dssr_stage_counts[1], 1)

    def test_directional_labels_and_insertion_codes(self):
        a, b = residue("AA", "-2"), residue("AA", "10A")
        self.assertEqual(ingest([(b, a, "cHW")]), {(a, b): {"cWH"}})
        self.assertEqual(flip_label("cWH|tHS"), "cHW|tSH")
        self.assertEqual(residue("AA", "10", "A"), b)

    def test_disagreement_is_not_correct_classification(self):
        metrics = self.c.metrics()
        self.assertEqual(metrics["cHS"]["matched"], 0)
        self.assertEqual(metrics["cHS"]["total_python"], 1)
        self.assertEqual(metrics["cWH"]["total_dssr"], 1)

    def test_mismatched_archives_excluded_from_both_stages(self):
        self.assertEqual(self.c.cohort, {"P1"})
        self.assertEqual(self.c.selected.keys(), self.c.identified.keys())
        self.assertLessEqual(self.c.selected["P1"].keys(), self.c.identified["P1"].keys())

    def test_unclassified_and_invalid_candidates_do_not_enter_counts(self):
        counts = self.c.disposition_counts()
        self.assertEqual(counts.sum(), 5)
        self.assertEqual(counts["selected"], 2)
        self.assertEqual(counts["ineligible"], 1)


if __name__ == "__main__":
    unittest.main()
