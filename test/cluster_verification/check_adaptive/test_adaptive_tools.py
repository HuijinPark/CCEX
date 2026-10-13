#!/usr/bin/env python3
"""Tests for the verification tools themselves.

Run: python -m unittest discover -s test/cluster_verification/check_adaptive -p 'test_*.py'
A comparison that could not be made must never look like "no difference", and only
elapsed measurements may be dropped from a comparison.
"""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import compare_summaries as cs


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    return path


class CompareSummaries(unittest.TestCase):
    def test_missing_file_is_an_error(self):
        with tempfile.TemporaryDirectory() as root:
            a = write(Path(root)/'a.json', {'x': 1}); missing = Path(root)/'missing.json'
            with self.assertRaises(cs.LoadError):
                cs.compare(a, missing)
            done = subprocess.run([sys.executable, str(HERE/'compare_summaries.py'), str(a), str(missing)],
                                  capture_output=True, text=True)
            self.assertEqual(done.returncode, 2)
            self.assertEqual(done.stdout, '')

    def test_invalid_json_is_an_error(self):
        with tempfile.TemporaryDirectory() as root:
            a = write(Path(root)/'a.json', {'x': 1})
            b = Path(root)/'b.json'; b.write_text('{')
            with self.assertRaises(cs.LoadError):
                cs.compare(a, b)

    def test_only_elapsed_measurements_are_dropped(self):
        settings = {'max_total_seconds': 30, 'seconds_limit': 120, 'raw_propagation_started': False,
                    'full_time_window': True}
        with tempfile.TemporaryDirectory() as root:
            a = write(Path(root)/'a.json', {'x': 1, 'wall_s': 3, 'search_wall_s': 1, 'search_s': .5, **settings})
            b = write(Path(root)/'b.json', {'x': 1, 'wall_s': 9, 'search_wall_s': 2, 'search_s': .9, **settings})
            self.assertEqual(cs.compare(a, b, drop=True), [])
            self.assertEqual(cs.compare(a, b), ['/search_s', '/search_wall_s', '/wall_s'])
            c = write(Path(root)/'c.json', {'x': 1, 'wall_s': 3, 'search_wall_s': 1, 'search_s': .5,
                                            **{**settings, 'max_total_seconds': 60,
                                               'raw_propagation_started': True}})
            self.assertEqual(cs.compare(a, c, drop=True), ['/max_total_seconds', '/raw_propagation_started'])


    def test_signed_zero_differs_and_nan_matches(self):
        with tempfile.TemporaryDirectory() as root:
            a = write(Path(root)/'a.json', {'x': 0.0, 'y': float('nan')})
            b = write(Path(root)/'b.json', {'x': -0.0, 'y': float('nan')})
            self.assertEqual(cs.compare(a, b), ['/x'])


class RegressVerdict(unittest.TestCase):
    """regress.py's comparison with golden (A46): decisions exact on every build; with MKL, floats within 1e-5 of the
    golden value, without it within 1e-9."""
    golden = {'returncode': 0, 'summaries': {'summary.json': {'rounds': [
        {'accepted': {'patch': [3, 7], 'policy': 'merge', 'local_after': 0.11628271974642508},
         'whole_bath_error': {'scaled_error': 0.11628271974641156}, 'stopped': False}]}}}

    def verdict(self, run, mkl):
        import regress
        return list(regress.differences(self.golden, run, mkl=mkl))

    def changed(self, path, value):
        run = json.loads(json.dumps(self.golden))
        node = run
        for k in path[:-1]:
            node = node[k]
        node[path[-1]] = value
        return run

    def test_a_different_decision_fails_with_mkl(self):
        base = ['summaries', 'summary.json', 'rounds', 0]
        for path, value in ((base + ['accepted', 'patch'], [3, 8]), (base + ['accepted', 'policy'], 'regroup'),
                            (base + ['stopped'], True), (['returncode'], 1)):
            self.assertTrue(self.verdict(self.changed(path, value), mkl=True), path)

    def test_floats_with_mkl(self):
        path = ['summaries', 'summary.json', 'rounds', 0, 'whole_bath_error', 'scaled_error']
        g = self.golden['summaries']['summary.json']['rounds'][0]['whole_bath_error']['scaled_error']
        self.assertEqual(self.verdict(self.changed(path, g*(1+4.4e-6)), mkl=True), [])     # the WSL MKL build's worst
        self.assertTrue(self.verdict(self.changed(path, g*(1+1e-4)), mkl=True))
        self.assertTrue(self.verdict(self.changed(path, g*(1+4.4e-6)), mkl=False))           # without MKL: 1e-9

    def test_zero_golden_is_absolute(self):
        import regress
        self.assertTrue(regress.close(0.0, 1e-13, True))
        self.assertFalse(regress.close(0.0, 1e-11, True))
        self.assertTrue(regress.close(1.0026641665932502e-15, 8.88185195937854e-16, True))   # round-off sized
        self.assertFalse(regress.close(1e-9, 3e-9, True))

    def test_mkl_detected_from_the_run_logs(self):
        import regress
        with tempfile.TemporaryDirectory() as d:
            out = Path(d)
            write(out/'0001-run'/'run.log', 'x').write_text('               MKL for large clusters : off\n')
            self.assertFalse(regress.built_with_mkl(out))
            write(out/'0002-probe'/'run.log', 'x').write_text('               MKL for large clusters : on (zgemm from dim 21)\n')
            self.assertTrue(regress.built_with_mkl(out))


class FaultLauncher(unittest.TestCase):
    def test_fixed_profile_rewrites_the_timings(self):
        with tempfile.TemporaryDirectory() as root:
            fake = Path(root)/'fake-mpiexec'
            fake.write_text('#!/bin/sh\necho "[profile] rank=0 method=cce backend=cpu dim=4 clusters=5'
                            ' build_s=0.5 setup_s=0.25 evolution_s=0.125"\necho done\n')
            fake.chmod(0o755)
            step = Path(root)/'0001-probe'; step.mkdir()
            (step/'input.json').write_text('{}')
            done = subprocess.run([sys.executable, str(HERE/'fault_mpiexec.py'), '-n', '1', 'binary',
                                   '-f', str(step/'input.json')], capture_output=True, text=True,
                                  env=dict(os.environ, CCEX_PROFILE='flat', CCEX_REAL_MPIEXEC=str(fake)))
            self.assertEqual(done.returncode, 0, done.stderr)
            fields = dict(part.split('=', 1) for part in done.stdout.splitlines()[0].split()[1:])
            self.assertEqual(float(fields['build_s']), 0.)
            self.assertEqual(float(fields['setup_s']), 0.)
            self.assertEqual(float(fields['evolution_s']), 5*1e-3)
            self.assertIn('done', done.stdout)

    def test_hanging_ranks_record_their_start(self):
        with tempfile.TemporaryDirectory() as root:
            fake = Path(root)/'fake-mpiexec'      # stands in for mpirun: runs the program once
            fake.write_text('#!/bin/sh\nshift 2\nexec "$@"\n')
            fake.chmod(0o755)
            step = Path(root)/'0005-probe'; step.mkdir()
            (step/'input.json').write_text('{}')
            log = Path(root)/'fault.jsonl'
            # The ranks sleep until killed, exactly as the timeout cases need.
            proc = subprocess.Popen([sys.executable, str(HERE/'fault_mpiexec.py'), '-n', '1', 'binary',
                                     '-f', str(step/'input.json')], start_new_session=True,
                                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                    env=dict(os.environ, CCEX_FAULT='0005-probe:hang-mpi',
                                             CCEX_REAL_MPIEXEC=str(fake), CCEX_FAULT_LOG=str(log)))
            try:
                deadline = time.monotonic()+10
                while time.monotonic() < deadline and not sorted(step.glob('hang-started.*')):
                    time.sleep(.05)
            finally:
                os.killpg(proc.pid, signal.SIGKILL); proc.wait()
            record = json.loads(log.read_text().splitlines()[0])
            self.assertEqual(record['action'], 'hang-mpi')
            self.assertTrue(record['marker'].startswith('ccex-hang-'))
            self.assertTrue(sorted(step.glob('hang-started.*')))


if __name__ == '__main__':
    unittest.main()
