"""Local, tiny-process tests for input preservation and the server resource guard."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

from run_orders import INPUTS, monitored_run, prepare_config


class ServerDriverTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    @staticmethod
    def resources():
        return dict(host_available_bytes=16 * 1024**3, host_total_bytes=32 * 1024**3,
                    gpu_used_mib=100, gpu_total_mib=96000)

    def test_original_physics_preserved(self):
        source = INPUTS / 'Diamond_NV_P1_10ppm_gcce/inputfiles/ccein_nstate1.json'
        original = json.loads(source.read_text())
        generated = prepare_config(source, self.root, 5, 'auto', 4)
        allowed = {'order', 'backend', 'profile', 'savemode', 'outfile', 'gpu_options',
                   'bathfile', 'qubitfile', 'gyrofile', 'avaaxfile', 'statefile', 'exstatefile'}
        for key, value in original.items():
            if key not in allowed:
                self.assertEqual(generated[key], value, key)
        self.assertEqual(generated['order'], 5)
        self.assertEqual(len(generated['bathfile']), len(original['bathfile']))
        for new, old in zip(generated['bathfile'], original['bathfile']):
            self.assertEqual(Path(new), (source.parent.parent / old).resolve())

    def test_success(self):
        result = monitored_run([sys.executable, '-c', 'print("finished")'], self.root,
                               5, sampler=self.resources, interval=0.01)
        self.assertIn('finished', (self.root / 'run.log').read_text())
        self.assertGreater(result['wall_s'], 0)

    def test_nonzero_exit(self):
        with self.assertRaisesRegex(RuntimeError, 'status 7'):
            monitored_run([sys.executable, '-c', 'raise SystemExit(7)'], self.root,
                          5, sampler=self.resources, interval=0.01)

    def test_memory_guard_stops_child(self):
        def low_memory():
            return dict(self.resources(), host_available_bytes=1)
        with self.assertRaises(MemoryError):
            monitored_run([sys.executable, '-c', 'import time; time.sleep(60)'], self.root,
                          5, sampler=low_memory, interval=0.01)

    def test_timeout_stops_child(self):
        with self.assertRaises(TimeoutError):
            monitored_run([sys.executable, '-c', 'import time; time.sleep(60)'], self.root,
                          0.05, sampler=self.resources, interval=0.01)


if __name__ == '__main__':
    unittest.main()

