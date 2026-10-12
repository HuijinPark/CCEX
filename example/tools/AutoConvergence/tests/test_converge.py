"""Regression cases of the review 2026-10-07 ((internal review)): A1, A2 (runner.converge) and A3
(window.fit_window), on stub curves: an area per value, rel = |A_lo/A_hi - 1|, pass within 2 %."""
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from autoconv import window                      # noqa: E402
from autoconv.runner import converge, unjudged_reason   # noqa: E402


def judge(a, b):
    if a is None or b is None:
        return dict(verdict="fit-failure", rel=None)
    rel = abs(a / b - 1)
    return dict(verdict="pass" if rel <= 0.02 else "fail", rel=rel)


class TestConverge(unittest.TestCase):
    def run_case(self, areas, values, default):
        return converge(values, lambda v: areas.get(v, default), lambda v: None, 10, False, judge, 0.02)

    def test_a_larger_failing_value_blocks_the_trend_choice(self):      # A1
        r = self.run_case({1: 97, 2: 99, 3: 100, 3.75: 97.1}, [1, 2, 3], 97.1)
        self.assertEqual(r["estimate"]["status"], "accepted")            # the trend still proposes 2 ...
        self.assertNotEqual(r["chosen"], 2)                              # ... but 3 fails against the larger value
        self.assertEqual(r["chosen"], 3.75)
        self.assertEqual([x["verdict"] for x in r["vs_reference"] if x["value"] == 3], ["fail"])

    def test_a_failed_fit_gives_no_trend_instead_of_a_crash(self):     # A2
        r = self.run_case({1: None, 2: 97, 3: 99, 4: 100}, [1, 2, 3, 4], 100.1)
        self.assertTrue(r["estimate"]["status"].startswith("no trend"))
        self.assertNotEqual(r["chosen"], 1)


class TestUnjudgedReason(unittest.TestCase):      # 2026-10-12, job 17163: a fit failure is not "make larger baths"
    def test_last_comparison_a_fit_failure(self):
        rec = dict(checks=[dict(lo=1, hi=2, verdict="fail", rel=0.5), dict(lo=2, hi=3, verdict="fit-failure", rel=None)])
        why = unjudged_reason("rbath", rec, "3_rbath.png")
        self.assertIn("could not be judged", why)
        self.assertIn("2 vs 3 (fit-failure)", why)
        self.assertNotIn("1 vs 2", why)                  # a judged comparison is not listed
        self.assertNotIn("larger baths", why)

    def test_last_comparison_not_finite(self):
        rec = dict(checks=[dict(lo=4, hi=5, verdict="nonfinite", rel=None)])
        self.assertIn("4 vs 5 (nonfinite)", unjudged_reason("rdip", rec, "2_rdip.png"))

    def test_judged_last_comparison_keeps_the_usual_reason(self):
        rec = dict(checks=[dict(lo=1, hi=2, verdict="fit-failure", rel=None), dict(lo=2, hi=3, verdict="fail", rel=0.1)])
        self.assertIsNone(unjudged_reason("rbath", rec, "3_rbath.png"))

    def test_nothing_compared(self):
        self.assertIsNone(unjudged_reason("rbath", None, "3_rbath.png"))
        self.assertIsNone(unjudged_reason("rbath", dict(checks=[]), "3_rbath.png"))


class TestPlateau(unittest.TestCase):
    def test_a_curve_near_its_start_is_not_a_plateau(self):            # A3
        t = np.linspace(0, 1, 101)
        data = [np.exp(-t / 100)[None, :].astype(complex) for _ in range(3)]
        self.assertEqual(window.fit_window(t, data, margin=1.5)["status"], "needs-extension")

    def test_a_decayed_curve_that_stops_is_a_plateau(self):
        t = np.linspace(0, 10, 201)
        data = [(0.3 + 0.7 * np.exp(-t))[None, :].astype(complex) for _ in range(3)]
        self.assertEqual(window.fit_window(t, data, margin=1.5)["status"], "plateau")


if __name__ == "__main__":
    unittest.main()
