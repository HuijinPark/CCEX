"""Tests of make_state_files.py (example/tools/StateFiles).

    cd example/tools/StateFiles && python3 -m unittest discover -s tests -v

The last test runs bin/main.out (the repository's build) on the files the tool writes and checks that CCEX reads every one
of them; it is skipped without a build, or with a build older than configuration.s_range (A94).
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import make_state_files as msf  # noqa: E402

ROOT = HERE.parents[3]
P1 = ROOT / "test/code_verification/bath_data/single_data"


def bath(path, names):
    Path(path).write_text(f"{len(names) + 1}\t0 0 0\n" + "".join(f"{i} {2 * i} {3 * i} {n}\n" for i, n in enumerate(names)))


class TestNames(unittest.TestCase):
    def test_k_and_s(self):
        self.assertEqual(msf.name_for("cfg{k}_state{s}", 3, 7), "cfg3_state7")
        self.assertEqual(msf.name_for("cfg{k}_state_", 3, 7), "cfg3_state_7")          # no {s}: the number appended
        self.assertEqual(msf.name_for("a{{s}}_{s}", None, 2), "a{s}_2")
        self.assertEqual(msf.name_for("a{{s}}", None, 2), "a{s}2")
        self.assertEqual(msf.name_for("a{{k}}{k}", 4), "a{k}4")
        self.assertEqual(msf.name_for("a{k}", None, None), "a{k}")                      # no k_range: {k} kept

    def test_comments_are_read_as_ccex_reads_them(self):
        text = '{"a": "x//y#z!", // c\n "b": 1 # c\n, "c": 2 ! c\n}'
        self.assertEqual(json.loads(msf.strip_comments(text)), {"a": "x//y#z!", "b": 1, "c": 2})


class TestFiles(unittest.TestCase):
    def setUp(self):
        self.d = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.d)

    def write_input(self, conf, species=None, name="in.jsonc"):
        inp = {"configuration": conf}
        if species is not None:
            inp["Defect"] = {"species": species}
        (self.d / name).write_text("// a test input\n" + json.dumps(inp, indent=1))
        return self.d / name

    def test_every_file_in_ccex_formats(self):
        (self.d / "data").mkdir()
        for k in (1, 2):
            bath(self.d / f"data/bath_{k}", ["13C", "P1", "13C", "P1", "14N"])
        (self.d / "data/gyro").write_text("P1 0.5 -17608.597050\n")
        species = [{"name": "p1", "n_axis_states": 4, "n_onsite_spins": 2, "types": ["14N", "13C"], "spins": [1.0, 0.5]}]
        inp = self.write_input({"path": "data", "gyrofile": "gyro", "k_range": [1, 2], "s_range": [2, 3], "seed": 5,
                                "bathfile": ["bath_{k}"], "bath_state_file": ["st/c{k}_s{s}"],
                                "defect_axis_file": ["ax_{k}"], "defect_state_file": ["ex_{k}_"]}, species)
        written = msf.make(inp)
        self.assertEqual(len(written), 2 * 2 + 2 + 2 * 2)
        for k in (1, 2):
            for s in (2, 3):
                st = (self.d / f"data/st/c{k}_s{s}").read_text().split("\n")
                self.assertEqual(st[0], "6")                                   # the bath file's first line
                m = [float(x) for x in st[1:6]]
                u = np.random.default_rng([5, k, s]).random(5)                 # the rule
                S = np.array([0.5, 0.5, 0.5, 0.5, 1.0])
                self.assertEqual(m, list(S - np.floor(u * (2 * S + 1))))
                ex = [l.split() for l in (self.d / f"data/ex_{k}_{s}").read_text().splitlines()]
                self.assertEqual([len(r) for r in ex], [1, 1, 2, 1, 2, 1])     # first line, then 13C, P1, 13C, P1, 14N
                self.assertTrue(all(float(r[0]) in (1, 0, -1) and float(r[1]) in (0.5, -0.5) for r in (ex[2], ex[4])))
                self.assertTrue(ex[1] == ex[3] == ex[5] == ["0"])
            ax = [int(x) for x in (self.d / f"data/ax_{k}").read_text().split()]
            self.assertEqual(ax[0], 6)
            self.assertTrue(ax[1] == ax[3] == ax[5] == 0 and 1 <= ax[2] <= 4 and 1 <= ax[4] <= 4)
        self.assertNotEqual((self.d / "data/st/c1_s2").read_text(), (self.d / "data/st/c1_s3").read_text())
        with self.assertRaises(msf.Refused):                                  # the files exist
            msf.make(inp)
        before = (self.d / "data/st/c1_s2").read_text()
        msf.make(inp, force=True)
        self.assertEqual((self.d / "data/st/c1_s2").read_text(), before)        # reproducible

    def test_the_runner_makes_the_same_files(self):                           # one rule (A94)
        sys.path.insert(0, str(ROOT / "example/tools/AutoConvergence"))
        try:
            from autoconv.runner import write_states, write_avaax
        finally:
            sys.path.pop(0)
        bath(self.d / "b_7", ["13C", "P1", "13C"])
        (self.d / "gyro").write_text("P1 0.5 -17608.6\n")
        inp = self.write_input({"gyrofile": "gyro", "k_range": [7, 7], "s_range": [1, 2], "seed": 3, "bathfile": ["b_{k}"],
                                "bath_state_file": ["s{k}_"], "defect_axis_file": ["a{k}"]},
                               [{"name": "P1", "n_axis_states": 3}])
        msf.make(inp)
        write_states(self.d / "r", self.d / "b_7", 7, 2, 3, {"13C": 0.5, "P1": 0.5})
        write_avaax(self.d / "ra", self.d / "b_7", 7, 3, {"p1": 3})
        for s in (1, 2):
            self.assertEqual((self.d / f"s7_{s}").read_text(), (self.d / f"r{s}").read_text())
        self.assertEqual((self.d / "a7").read_text(), (self.d / "ra").read_text())

    def test_seeds_and_refusals(self):
        bath(self.d / "b1", ["13C"]); bath(self.d / "b2", ["13C"])
        inp = self.write_input({"k_range": [1, 2], "s_range": [1, 1], "seed": [11, 12], "bathfile": ["b{k}"],
                                "bath_state_file": ["s{k}_"]})
        self.assertEqual([x[5] for x in msf.plan(msf.read_input(inp), self.d)], [11, 12])
        plan = msf.plan(msf.read_input(self.write_input({"bathfile": ["b1", "b2"], "s_range": [1, 1],
                                                         "bath_state_file": ["", "s_"]})), self.d)
        self.assertEqual([(str(p.name), k, x) for _, p, _, k, _, x in plan], [("s_1", 0, 1)])   # "" skipped, k 0, seed 1
        for conf, msg in (({"bathfile": ["b1"], "nstate": 2}, "s_range"),
                          ({"bathfile": ["b1"], "s_range": [0, 2]}, "s_range"),
                          ({"bathfile": ["b1"], "s_range": 2}, "s_range"),
                          ({"bathfile": ["b{k}"], "k_range": [1, 2], "seed": [1]}, "one integer per configuration"),
                          ({"bathfile": ["b1"], "bath_state_file": ["a", "b"], "s_range": [1, 1]}, "one name per bath file")):
            with self.subTest(conf=conf), self.assertRaisesRegex(msf.Refused, msg):
                msf.plan(msf.read_input(self.write_input(conf)), self.d)
        bath(self.d / "bx", ["Xx"])
        with self.assertRaisesRegex(msf.Refused, "neither in the gyrofile"):
            msf.make(self.write_input({"bathfile": ["bx"], "s_range": [1, 1], "bath_state_file": ["q"]}, name="x.json"))

    def test_command_line(self):
        bath(self.d / "b1", ["13C", "13C"])
        inp = self.write_input({"bathfile": ["b1"], "s_range": [1, 3], "bath_state_file": ["st{s}"]})
        run = lambda *a: subprocess.run([sys.executable, str(HERE.parent / "make_state_files.py"), *a], capture_output=True, text=True)
        p = run(str(inp))
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("3 files written", p.stdout)
        self.assertEqual(run(str(inp)).returncode, 2)
        self.assertEqual(run(str(inp), "--force").returncode, 0)


@unittest.skipUnless((ROOT / "bin/main.out").is_file() and (P1 / "bath_DiaP1_1ppm_1").is_file(), "no bin/main.out or no P1 bath")
class TestCcexReadsThem(unittest.TestCase):
    def test_a_run_reads_every_file(self):
        d = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, d)
        shutil.copy(P1 / "bath_DiaP1_1ppm_1", d / "bath_1")
        shutil.copy(P1 / "DiaP1_gyro", d / "gyro")
        inp = {"method": "cCCE", "order": 2, "bfield": 500, "rbath": 400, "rdip": 400, "deltat": 0.01, "nstep": 3,
               "pulse": {"count": 1},
               "Qubit": {"qubit": [{"file": str(P1 / "bath_DiaP1_1ppm_Defect"), "spin": 1.0}]},
               "configuration": {"gyrofile": "gyro", "bathfile": ["bath_{k}"], "k_range": [1, 1], "s_range": [1, 2], "seed": 4,
                                 "bath_state_file": ["state_{k}_{s}"], "defect_axis_file": ["axis_{k}"],
                                 "defect_state_file": ["onsite_{k}_"]},
               "Defect": {"species": [{"name": "P1", "n_onsite_spins": 2, "n_axis_states": 4, "types": ["14N", "13C"],
                                       "spins": [1.0, 0.5], "gyros": [1.9338, 6.728284]}]}}   # two on-site spins: one line holds two values
        (d / "in.json").write_text(json.dumps(inp, indent=1))
        msf.make(d / "in.json")
        p = subprocess.run([str(ROOT / "bin/main.out"), "-f", "in.json"], cwd=d, capture_output=True, text=True,
                           env=dict(os.environ, OMP_NUM_THREADS="1"))
        log = p.stdout + p.stderr
        if p.returncode != 0 and any(w in log for w in ("s_range", "Qubit.file", "nstate", "xyz is not found")):
            self.skipTest("bin/main.out is older than configuration.s_range and Qubit.qubit[].file (A94)")
        self.assertEqual(p.returncode, 0, log[-3000:])
        for word in ("Read StateFile : state_1_1", "Read StateFile : state_1_2", "Read AvaaxFile : axis_1",
                     "Read ExstateFile : onsite_1_1", "Read ExstateFile : onsite_1_2"):
            self.assertIn(word, log)
        self.assertNotIn("randomly generated", log)
        self.assertTrue((d / "in_k1_state2_wiDiv").is_file() or any(d.glob("*state2_wiDiv")), sorted(x.name for x in d.iterdir()))


if __name__ == "__main__":
    unittest.main()
