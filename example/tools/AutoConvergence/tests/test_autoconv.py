"""Tests of the auto_convergence setup driver (time window and sampling steps).

    cd example/tools/AutoConvergence && python -m unittest discover -s tests -v

Synthetic cases check the estimators against known answers; the last test reads the existing CR0 runs of the
auto-truncation study (~/tmp/ccex-auto-trunc/0c-res/A) when they are present and is skipped otherwise.
"""
import contextlib
import io
import json
import math
import re
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from autoconv import config, curves, sampling, window  # noqa: E402
from autoconv.__main__ import _finite, analyze, main  # noqa: E402


def write_run(d, times, states, prefix="out"):
    d.mkdir(parents=True, exist_ok=True)
    for s, v in enumerate(states, 1):
        lines = [f"{t:15g}\t{x.real:+15.10f}{x.imag:+.10f}j\t" for t, x in zip(times, v)]
        (d / f"{prefix}_state{s}_wiDiv").write_text("\n".join(lines) + "\n")
        (d / f"{prefix}_state{s}_noDiv").write_text("\n".join(lines) + "\n")   # must not be read


class TestConfig(unittest.TestCase):
    def test_defaults_and_budget(self):
        s = config.parse_auto_convergence({})
        self.assertEqual(s.tolerance, 0.01)
        self.assertAlmostEqual(s.rdip_tolerance, 0.01 / 3 * 0.5)
        self.assertAlmostEqual(s.rbath_tolerance, 0.01 / 3 * 0.5)
        self.assertAlmostEqual(s.sampling_tolerance, 0.01 * 2 / 3)
        self.assertLessEqual(s.rdip_tolerance + s.rbath_tolerance + s.sampling_tolerance, s.tolerance * (1 + 1e-12))
        self.assertEqual(s.steps, {"time_window": True, "rdip": False, "rbath": False, "sampling": True})

    def test_refusals(self):
        bad = [{"tau": 1e-3}, {"steps": {"rdip": True}}, {"rdip_start": 200}, {"error_norm": "l2"},
               {"tolerance": -1}, {"decay_threshold": 1.5}, {"window_margin": 0.5}, {"steps": {"time_window": 1}},
               {"tolerance": 0.01, "sampling_tolerance": 0.02}]
        for b in bad:
            with self.subTest(b=b), self.assertRaises(config.ConfigError):
                config.parse_auto_convergence(b)

    def test_base_input_consistency(self):
        good = {"rdip": 300, "rbath": 288, "deltat": 1e-4, "nstep": 151}
        config.parse_auto_convergence({}, good)
        for drop in ("rdip", "rbath", "deltat"):
            with self.subTest(drop=drop), self.assertRaises(config.ConfigError):
                config.parse_auto_convergence({}, {k: v for k, v in good.items() if k != drop})

    def test_load_config(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "ccein.json").write_text('{"rdip": 300, "rbath": 288, # comment\n "deltat": 1e-4, "nstep": 151}')
            drv = {"base_input": "ccein.json", "auto_convergence": {"tolerance": 0.02},
                   "bath_source": {"files": "b*"}, "state_source": {"files": "s{k}"}}
            (tmp / "drv.json").write_text(json.dumps(drv))
            s = config.load_config(tmp / "drv.json")
            self.assertEqual(s.tolerance, 0.02)
            self.assertEqual(s.base_input["rbath"], 288)
            for broken in ({**drv, "bath_source": {}}, {**drv, "state_source": {}}, {**drv, "extra": 1},
                           {**drv, "bath_source": {"files": "b*", "generator": "g.py"}}):
                (tmp / "drv.json").write_text(json.dumps(broken))
                with self.subTest(broken=broken), self.assertRaises(config.ConfigError):
                    config.load_config(tmp / "drv.json")


class TestCurves(unittest.TestCase):
    def test_read_run_and_refusals(self):
        t = np.linspace(0, 1, 11)
        states = [np.exp(-t) * (1 + 0.1j * s) for s in range(3)]
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "r1"
            write_run(d, t, states)
            _, x = curves.read_run(d)
            self.assertEqual(x.shape, (3, 11))
            self.assertTrue(np.allclose(x, np.array(states), atol=1e-9))
            (d / "out_state2_wiDiv").unlink()
            with self.assertRaises(ValueError):
                curves.read_run(d)
            d2 = Path(tmp) / "r2"
            write_run(d2, t, states)
            write_run(Path(tmp) / "r3", t * 2, states)
            with self.assertRaises(ValueError):
                curves.read_runs([d2, Path(tmp) / "r3"])


class TestSampling(unittest.TestCase):
    def synthetic(self, c=400, s=10, sb=0.05, sw=0.2, nt=40, seed=1):
        rng = np.random.default_rng(seed)
        cplx = lambda *shape: (rng.normal(size=shape) + 1j * rng.normal(size=shape)) / math.sqrt(2)
        base = np.exp(-np.linspace(0, 3, nt))
        cfg = base + sb * cplx(c, nt)
        return [cfg[i] + sw * cplx(s, nt) for i in range(c)]

    def test_variance_components_recover_truth(self):
        comp = sampling.variance_components(self.synthetic())
        self.assertAlmostEqual(float(np.sqrt(comp["W"]).mean()), 0.2, delta=0.01)
        self.assertAlmostEqual(float(np.sqrt(comp["B"]).mean()), 0.05, delta=0.008)

    def test_configurations_needed_is_minimal(self):
        comp = sampling.variance_components(self.synthetic())
        for s in (1, 5, 30):
            c = sampling.configurations_needed(comp, s, 0.003)
            self.assertLessEqual(sampling.sampling_rms_proxy(comp, c, s), 0.003 * (1 + 1e-12))
            if c > 2:
                self.assertGreater(sampling.sampling_rms_proxy(comp, c - 1, s), 0.003)

    def test_allocation_prefers_more_configurations_when_states_cost_as_much(self):
        comp = sampling.variance_components(self.synthetic())
        cheap_cfg = sampling.allocate(comp, 0.003, cost_config=0.0, cost_state=1.0)["best"]
        dear_cfg = sampling.allocate(comp, 0.003, cost_config=50.0, cost_state=1.0)["best"]
        self.assertLessEqual(cheap_cfg["S"], dear_cfg["S"])

    def test_bootstrap_shrinks_with_configurations(self):
        few = sampling.bootstrap_norm_error(self.synthetic(c=10), n_boot=500)
        many = sampling.bootstrap_norm_error(self.synthetic(c=200), n_boot=500)
        self.assertLess(many["mean"], few["mean"])

    def test_needs_two_states_and_configurations(self):
        with self.assertRaises(ValueError):
            sampling.variance_components([np.ones((1, 5))] * 3)
        with self.assertRaises(ValueError):
            sampling.variance_components([np.ones((3, 5))])


class TestWindow(unittest.TestCase):
    t = np.linspace(0, 10, 201)

    def configs(self, curve, c=8, s=5, noise=1e-4, seed=2):
        rng = np.random.default_rng(seed)
        return [curve + noise * rng.normal(size=(s, len(curve))) for _ in range(c)]

    def test_revival_after_first_crossing_is_not_the_decay(self):
        curve = np.exp(-self.t) + 0.05 * np.exp(-((self.t - 6.0) / 0.3) ** 2)   # dips below 0.01 near t=4.6, revives at 6
        rep = window.choose_window(self.t, self.configs(curve), threshold=0.01, margin=1.0)
        self.assertLess(rep["first_crossing"], 5.0)
        self.assertGreater(rep["max_after_first_crossing"], 0.04)
        self.assertGreater(rep["decay_time"], 6.3)            # the revival exceeds 0.01 until ~6.38
        self.assertEqual(rep["status"], "resolved")

    def test_status_needs_extension_and_unresolved(self):
        rep = window.choose_window(self.t, self.configs(np.exp(-self.t)), threshold=0.01, margin=1.5)
        self.assertAlmostEqual(rep["decay_time"], math.log(100), delta=0.1)
        self.assertEqual(rep["status"], "resolved")                       # 1.5 * 4.6 = 6.9 < 10
        rep = window.choose_window(self.t, self.configs(np.exp(-self.t / 3)), threshold=0.01, margin=1.5)
        self.assertEqual(rep["status"], "window-unresolved")             # exp(-10/3) = 0.036 never below 0.01
        rep = window.choose_window(self.t, self.configs(np.exp(-self.t / 1.9)), threshold=0.01, margin=1.5)
        self.assertEqual(rep["status"], "needs-extension")               # T_d ~ 8.7, needs 13 > 10

    def test_uses_complex_mean_not_mean_modulus(self):
        # states with random phases: |mean| decays fast although every |L| stays 1
        rng = np.random.default_rng(3)
        data = [np.exp(1j * rng.uniform(0, 2 * np.pi, (50, 1)) * self.t) for _ in range(20)]
        mu, _ = window.pooled_mean(data)
        self.assertLess(abs(mu[-1]), 0.2)

    def test_tail_does_not_dilute_the_sampling_error(self):              # 2026-10-07: area on [0, T_d]
        rng = np.random.default_rng(5)
        curves_ = [[np.exp(-6 * np.linspace(0, 2, 81)) * (1 + 0.1 * rng.normal()) + 0j for _ in range(4)] for _ in range(6)]
        s = config.parse_auto_convergence({})
        reps = []
        with tempfile.TemporaryDirectory() as tmp:
            for points in (41, 81):                      # the same curves on [0, 1] and on [0, 2] (a long near-zero tail added)
                runs = []
                for c, states in enumerate(curves_):
                    d = Path(tmp) / f"w{points}c{c}"
                    write_run(d, np.linspace(0, 2, 81)[:points], [x[:points] for x in states])
                    runs.append(str(d))
                reps.append(analyze(runs, s)["sampling"])
        self.assertIn("[0, T_d]", reps[0]["evaluated_on"])
        self.assertAlmostEqual(reps[0]["sampling_rms_proxy"], reps[1]["sampling_rms_proxy"], places=12)

    def test_grid_tolerance_key_is_gone(self):                           # 2026-10-07: nstep is the user's input
        with self.assertRaises(config.ConfigError):
            config.parse_auto_convergence({"grid_tolerance": 1e-4})


class TestCommandLine(unittest.TestCase):
    def test_analyze_end_to_end(self):
        t = np.linspace(0, 10, 101)
        rng = np.random.default_rng(4)
        with tempfile.TemporaryDirectory() as tmp:
            runs = []
            for c in range(6):
                d = Path(tmp) / f"cfg{c}"
                write_run(d, t, [np.exp(-t * (1 + 0.1 * rng.normal())) for _ in range(4)])
                runs.append(str(d))
            rep = analyze(runs, config.parse_auto_convergence({}), cost_config=2.0, cost_state=1.0)
            self.assertEqual(rep["configurations"], 6)
            self.assertIn(rep["time_window"]["status"], ("resolved", "needs-extension", "window-unresolved"))
            self.assertIn("allocation", rep["sampling"])
            out, err = io.StringIO(), io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                self.assertEqual(main(["analyze", *runs]), 0)
                self.assertEqual(main(["analyze", *runs, "--cost-config", "1"]), 2)
            self.assertEqual(json.loads(out.getvalue())["configurations"], 6)
            self.assertIn("both --cost-config and --cost-state", err.getvalue())


class TestReviewFixes(unittest.TestCase):
    """Regression cases of review 20261007-013640 (R2-R6)."""

    def test_unbalanced_states_use_one_sample(self):                    # R2
        t = np.linspace(0, 1, 3)
        a = np.array([[1, 0, 0], [1, 0, 0]], complex)
        b = np.array([[1, 0, 0], [1, 0, 0], [1, 1, 1]], complex)       # the third state changes the mean
        with tempfile.TemporaryDirectory() as tmp:
            write_run(Path(tmp) / "c1", t, a)
            write_run(Path(tmp) / "c2", t, b)
            rep = analyze([str(Path(tmp) / "c1"), str(Path(tmp) / "c2")], config.parse_auto_convergence({"steps": {"time_window": False}}, None))
        s = rep["sampling"]
        self.assertEqual((s["C"], s["S"], s["states_excluded"]), (2, 2, 1))
        self.assertEqual(s["bootstrap"]["quantile"], 0.0)                # bootstrap on the same 2 + 2 states as the proxy

    def test_bad_times_refused_non_finite_values_kept(self):              # R3, review 20261007d R8
        with tempfile.TemporaryDirectory() as tmp:
            for name, lines in (("inf_t", ["0 1+0j", "inf 0+0j", "inf 0+0j"]), ("order", ["0 1+0j", "2 0+0j", "1 0+0j"]),
                                ("empty", []), ("short", ["0", "1"])):
                f = Path(tmp) / name
                f.write_text("\n".join(lines) + "\n")
                with self.subTest(name=name), self.assertRaises(ValueError):
                    curves.read_curve(f)
            f = Path(tmp) / "split"
            f.write_text("0 1+0j\n1 -inf+inf j\n2 nan+0j\n")      # CCEX prints non-finite values with a space
            _, v = curves.read_curve(f)
            self.assertTrue(np.isinf(v[1].real) and np.isinf(v[1].imag) and np.isnan(v[2].real))

    def test_report_is_strict_json(self):                                 # R4
        comp = dict(W=np.ones(4), B=np.zeros(4), V=np.ones(4), negative_B=0, C=3, S=3)
        alloc = sampling.allocate(comp, 0.1, cost_config=0.0, cost_state=1.0)
        self.assertIsNone(alloc["S_star_continuous"])
        comp["B"] = np.ones(4)
        alloc = sampling.allocate(comp, 0.1, cost_config=0.0, cost_state=1.0)
        self.assertEqual(alloc["S_star_continuous"], 1.0)                 # a = 0 -> the smallest state count, not infinity
        json.dumps(alloc, allow_nan=False)

    def test_bang_and_hash_comments(self):                                # R5
        text = '{"rdip": 300, ! comment\n "rbath": 288, # other\n "path": "a#b!c", "deltat": 1e-4, "nstep": 151}'
        self.assertEqual(json.loads(config.strip_comments(text))["path"], "a#b!c")

    def test_small_time_scale_grids_compared(self):                       # R6
        self.assertFalse(curves.same_grid(np.array([0, 1e-9, 2e-9]), np.array([0, 2e-9, 4e-9])))
        self.assertTrue(curves.same_grid(np.array([0, 1.33333e-4]), np.array([0, 1.333333e-4])))


class TestReview20261007d(unittest.TestCase):
    """Regression cases of the plan review 20261007d (sample identity, numerical failures, allocation, extra states)."""

    def runs(self, tmp, n=3, states=3, points=5, extra=None):
        t = np.linspace(0, 1, points)
        rng = np.random.default_rng(7)
        out = []
        for c in range(n):
            k = states + (extra if extra and c == 0 else 0)
            d = Path(tmp) / f"cfg{c}"
            write_run(d, t, [np.exp(-t) * (1 + 0.05 * rng.normal()) for _ in range(k)])
            out.append(str(d))
        return out

    def test_same_configuration_twice_refused(self):                    # R5
        with tempfile.TemporaryDirectory() as tmp:
            r = self.runs(tmp)
            with self.assertRaises(ValueError):
                curves.read_runs([r[0], r[1], r[0]])
            import shutil
            shutil.copytree(r[0], Path(tmp) / "copy")
            with self.assertRaises(ValueError):
                curves.read_runs([r[0], str(Path(tmp) / "copy")])

    def test_incomplete_runs_refused(self):                              # R5
        with tempfile.TemporaryDirectory() as tmp:
            r = self.runs(tmp)
            curves.read_runs(r, expected_states=3, expected_points=5)
            for kw in (dict(expected_states=4), dict(expected_points=6)):
                with self.subTest(kw=kw), self.assertRaises(ValueError):
                    curves.read_runs(r, **kw)

    def test_non_finite_state_is_a_numerical_failure(self):              # R8
        with tempfile.TemporaryDirectory() as tmp:
            r = self.runs(tmp)
            f = Path(r[1]) / "out_state2_wiDiv"
            lines = f.read_text().splitlines()
            lines[3] = lines[3].split("\t")[0] + "\t-inf+inf j"
            f.write_text("\n".join(lines) + "\n")
            rep = analyze(r, config.parse_auto_convergence({}))
        self.assertEqual(rep["status"], "numerical-failure")
        self.assertEqual(rep["failed_configurations"], [r[1]])
        self.assertEqual(rep["configurations"], 3)                       # nothing dropped
        self.assertEqual(rep["validity"]["configurations"][1]["nonfinite_states"], [2])
        self.assertNotIn("sampling", rep)
        json.dumps(_finite(rep), allow_nan=False)

    def test_overflow_is_a_numerical_failure(self):                     # R8 (cfg64: finite but 1e286)
        with tempfile.TemporaryDirectory() as tmp:
            t = np.linspace(0, 1, 5)
            write_run(Path(tmp) / "a", t, [np.full(5, 1e286 + 0j), np.full(5, -1e286 + 0j)])
            write_run(Path(tmp) / "b", t, [np.full(5, 0.5 + 0j), np.full(5, 0.4 + 0j)])
            rep = analyze([str(Path(tmp) / "a"), str(Path(tmp) / "b")], config.parse_auto_convergence({"steps": {"time_window": False}}))
        self.assertEqual(rep["status"], "numerical-failure")
        self.assertEqual(rep["validity"]["configurations"][0]["U"], 2)

    def test_unphysical_states_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = self.runs(tmp)
            write_run(Path(r[2]), np.linspace(0, 1, 5), [np.full(5, 1.2 + 0j), np.full(5, 0.5 + 0j), np.full(5, 0.4 + 0j)])
            rep = analyze(r, config.parse_auto_convergence({}))
        self.assertEqual(rep["validity"]["configurations"][2]["unphysical_states"], [1])
        self.assertIn("warnings", rep)
        self.assertEqual(rep["total"]["status"], "unresolved")
        self.assertIn("rdip", rep["total"]["unmeasured"])

    def test_extra_states_used_for_the_within_check(self):              # R4
        with tempfile.TemporaryDirectory() as tmp:
            rep = analyze(self.runs(tmp, extra=7), config.parse_auto_convergence({"steps": {"time_window": False}}))
        s = rep["sampling"]
        self.assertEqual((s["S"], s["states_excluded"]), (3, 7))
        self.assertIn("within_sd_mean_all_states", s)

    def test_allocation_bounded_by_the_source(self):                    # R7
        comp = dict(W=np.full(4, 0.02), B=np.full(4, 0.0016), V=np.ones(4), negative_B=0, C=10, S=30)
        free = sampling.allocate(comp, 0.0067, 1.0, 1.0)
        self.assertGreater(free["best"]["C"], 70)
        self.assertNotIn(1, [r["S"] for r in free["table"]])            # the analysis needs 2 states per configuration
        capped = sampling.allocate(comp, 0.0067, 1.0, 1.0, max_configurations=70)
        self.assertTrue(capped["best"]["C"] <= 70)
        tight = sampling.allocate(comp, 0.0001, 1.0, 1.0, max_configurations=70)
        self.assertTrue(tight["status"].startswith("infeasible"))
        self.assertGreater(tight["best"]["shortfall"], 0)
        json.dumps(tight, allow_nan=False)

    def test_source_capacity_and_expected_shape(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            for k in range(4):
                (tmp / f"bath{k}").write_text("x")
            (tmp / "ccein.json").write_text('{"rdip": 300, "rbath": 288, "deltat": 1e-4, "nstep": 151, "configuration": {"nstate": 30}}')
            (tmp / "drv.json").write_text(json.dumps({"base_input": "ccein.json", "bath_source": {"files": "bath*"},
                                                      "state_source": {"files": "s{k}"}}))
            s = config.load_config(tmp / "drv.json")
            self.assertEqual(config.source_capacity(s), 4)
            self.assertEqual(config.expected_run_shape(s), (30, 151))


class TestReview20261007e(unittest.TestCase):
    """Regression cases of the plan re-review 20261007e (mixed pilot, manifest identity, grid, overflow, capacity, costs)."""

    def make(self, tmp, sizes, points=5, t=None, scale=None):
        t = np.linspace(0, 1, points) if t is None else t
        rng = np.random.default_rng(11)
        runs = []
        for c, k in enumerate(sizes):
            d = Path(tmp) / f"cfg{c}"
            write_run(d, t, [np.exp(-t) * (1 + 0.05 * rng.normal()) for _ in range(k)])
            runs.append(str(d))
        return runs

    def test_mixed_size_pilot_reads_with_per_run_counts(self):          # F1
        with tempfile.TemporaryDirectory() as tmp:
            runs = self.make(tmp, [10, 30, 10])
            _, data, _ = curves.read_runs(runs, expected_states=[10, 30, 10])
            self.assertEqual([x.shape[0] for x in data], [10, 30, 10])
            with self.assertRaises(ValueError):
                curves.read_runs(runs, expected_states=[10, 30, 30])      # a missing extension state is refused

    def test_manifest_ids_turn_identical_bytes_into_a_warning(self):   # F2
        import shutil
        with tempfile.TemporaryDirectory() as tmp:
            runs = self.make(tmp, [3, 3])
            shutil.copytree(runs[0], Path(tmp) / "copy")
            runs.append(str(Path(tmp) / "copy"))
            with self.assertRaises(ValueError):
                curves.read_runs(runs)
            _, _, warn = curves.read_runs(runs, manifest_ids=["a", "b", "c"])
            self.assertEqual(len(warn), 1)
            with self.assertRaises(ValueError):
                curves.read_runs(runs, manifest_ids=["a", "b", "a"])

    def test_wrong_grid_refused(self):                                  # F2
        with tempfile.TemporaryDirectory() as tmp:
            runs = self.make(tmp, [2, 2], points=3, t=np.array([5.0, 15.0, 25.0]))
            with self.assertRaises(ValueError):
                curves.read_runs(runs, expected_times=1.0 * np.arange(3))
            curves.read_runs(runs, expected_times=np.array([5.0, 15.0, 25.0]))

    def test_extra_state_overflow_is_a_numerical_failure(self):         # F4
        with tempfile.TemporaryDirectory() as tmp:
            t = np.linspace(0, 1, 3)
            write_run(Path(tmp) / "a", t, [np.full(3, 0.5 + 0j), np.full(3, 0.5 + 0j), np.full(3, 1e286 + 0j)])
            write_run(Path(tmp) / "b", t, [np.full(3, 0.4 + 0j), np.full(3, 0.4 + 0j)])
            rep = analyze([str(Path(tmp) / "a"), str(Path(tmp) / "b")], config.parse_auto_convergence({"steps": {"time_window": False}}))
        self.assertEqual((rep["status"], rep["failed_stage"]), ("numerical-failure", "sampling"))

    def test_window_overflow_with_sampling_off(self):                   # F4, F9
        with tempfile.TemporaryDirectory() as tmp:
            t = np.linspace(0, 1, 3)
            write_run(Path(tmp) / "a", t, [np.full(3, 1e286 + 0j), np.full(3, 1e286 + 0j)])
            write_run(Path(tmp) / "b", t, [np.full(3, -1e286 + 0j), np.full(3, -1e286 + 0j)])
            rep = analyze([str(Path(tmp) / "a"), str(Path(tmp) / "b")], config.parse_auto_convergence({"steps": {"sampling": False}}))
        self.assertEqual((rep["status"], rep["failed_stage"]), ("numerical-failure", "time_window"))

    def test_sampling_off_listed_as_unmeasured(self):                   # F9
        with tempfile.TemporaryDirectory() as tmp:
            rep = analyze(self.make(tmp, [3, 3]), config.parse_auto_convergence({"steps": {"sampling": False}}))
        self.assertIn("sampling (fixed, not assessed)", rep["total"]["unmeasured"])

    def test_window_uses_the_balanced_sample(self):                     # F6
        t = np.linspace(0, 2, 3)
        with tempfile.TemporaryDirectory() as tmp:
            write_run(Path(tmp) / "a", t, [np.array([1, 0, 0], complex)] * 2)
            write_run(Path(tmp) / "b", t, [np.array([1, 0, 0], complex)] * 2 + [np.array([1, 1, 1], complex)])
            rep = analyze([str(Path(tmp) / "a"), str(Path(tmp) / "b")],
                          config.parse_auto_convergence({"window_method": "envelope"}))   # 3 points: nothing to fit
        self.assertEqual(rep["time_window"]["status"], "resolved")       # the extra [1,1,1] state does not move the decision

    def test_allocation_boundaries(self):                                # F3, F5
        comp = dict(W=np.full(4, 0.02), B=np.full(4, 0.0016), V=np.ones(4), negative_B=0, C=10, S=30)
        self.assertTrue(sampling.allocate(comp, 0.0067, 1.0, 1.0, max_configurations=1)["status"].startswith("source-limited"))
        self.assertTrue(sampling.allocate(comp, 0.0067, 1.0, 1.0, max_configurations=0)["status"].startswith("source-limited"))
        for bad in ((float("nan"), 1.0), (1.0, float("inf"))):
            with self.assertRaises(ValueError):
                sampling.allocate(comp, 0.0067, *bad)
        huge = sampling.allocate(comp, 1e-300, 1.0, 1.0)
        self.assertTrue(huge["status"].startswith("infeasible"))
        json.dumps(huge, allow_nan=False)

    def test_state_source_caps_the_state_count(self):                    # review 20261007f F4
        comp = dict(W=np.full(4, 0.04), B=np.zeros(4), V=np.ones(4), negative_B=0, C=20, S=10)
        free = sampling.allocate(comp, 0.0067, 4.0, 1.6, max_configurations=50)
        self.assertEqual(free["best"]["S"], 60)
        capped = sampling.allocate(comp, 0.0067, 4.0, 1.6, max_configurations=50, max_states=30)
        self.assertLessEqual(max(r["S"] for r in capped["table"]), 30)
        self.assertLessEqual(capped["best"]["S"], 30)

    def test_capacity_counts_existing_files_absolute_globs(self):       # F5
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "b1").write_text("x"); (tmp / "b2").write_text("x"); (tmp / "bdir").mkdir()
            s = config.parse_auto_convergence({})
            s.root = tmp
            s.bath_source = {"files": str(tmp / "b*")}
            self.assertEqual(config.source_capacity(s), 2)
            s.bath_source = {"files": ["b1", "./b1", "missing"]}
            self.assertEqual(config.source_capacity(s), 1)

    def test_cli_with_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            runs = self.make(tmp, [3, 5, 3])
            man = dict(configurations=[dict(id=f"c{i}", run=r, states=n) for i, (r, n) in enumerate(zip(runs, (3, 5, 3)))],
                       available_configurations=40, grid=dict(deltat=0.25, nstep=5))
            (Path(tmp) / "m.json").write_text(json.dumps(man))
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                self.assertEqual(main(["analyze", "--manifest", str(Path(tmp) / "m.json"), "--cost-config", "1", "--cost-state", "1"]), 0)
            rep = json.loads(out.getvalue())
            self.assertEqual(rep["provenance"], "manifest")
            self.assertEqual(rep["sampling"]["allocation"]["max_configurations"], 40)
            man["grid"]["deltat"] = 0.3
            (Path(tmp) / "m.json").write_text(json.dumps(man))
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(main(["analyze", "--manifest", str(Path(tmp) / "m.json")]), 2)


ROOT = Path.home() / "tmp/ccex-auto-trunc/0c-res/A"


@unittest.skipUnless(ROOT.is_dir(), "existing CR0 runs of the auto-truncation study not present")
class TestExistingCR0(unittest.TestCase):
    def test_matches_the_numbers_in_the_plan(self):
        runs = []
        for d in sorted(ROOT.glob("cfg*")):
            rungs = sorted((x for x in d.iterdir() if re.fullmatch(r"r[\d.]+", x.name)), key=lambda x: float(x.name[1:]))
            runs.append(rungs[-1])                     # the all-pairs rung of each configuration
        times, data, _ = curves.read_runs(runs)
        comp = sampling.variance_components(data)
        self.assertEqual((comp["C"], comp["S"]), (10, 30))
        self.assertAlmostEqual(float(np.sqrt(comp["W"]).mean()), 0.142, delta=0.002)
        self.assertAlmostEqual(float(np.sqrt(comp["B"]).mean()), 0.040, delta=0.002)
        self.assertEqual(comp["negative_B"], 28)
        self.assertAlmostEqual(sampling.sampling_rms_proxy(comp, 70, 30), 0.00579, delta=0.0001)
        rep = window.choose_window(times, data, 0.01, 1.5)
        self.assertAlmostEqual(rep["first_crossing"], 10.4e-3, delta=1e-4)   # ms; plan 0.6: pooled first crossing 10.4 us


if __name__ == "__main__":
    unittest.main()
