"""Tests of set_input_keys.py (A63).

    python3 -m unittest discover -s example/tools/InputKeys/tests
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

TOOL = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'set_input_keys.py')

INPUT = """{
    # a comment
    "method": "cCCE",   ! another
    "outfile": "run#1",
    "bathfile": ["../bath/bath_"],
    "Defect": {"species": []}
}
"""


class SetInputKeys(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.src = os.path.join(self.dir.name, 'in.json')
        self.out = os.path.join(self.dir.name, 'out.json')
        with open(self.src, 'w') as f:
            f.write(INPUT)

    def tearDown(self):
        self.dir.cleanup()

    def run_tool(self, *args, out=None):
        return subprocess.run([sys.executable, TOOL, self.src, out or self.out, *args], capture_output=True, text=True)

    def test_values_and_dotted_keys(self):
        p = self.run_tool('s_range=[1, 50]', 'bathfile=["../bath/bath_1"]', 'defect.axis_file=Random3', 'outfile=name_1',
                          'cluster_only=true')
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(self.out) as f:
            d = json.load(f)                 # plain JSON: the comments are gone
        self.assertEqual(d['s_range'], [1, 50])
        self.assertEqual(d['bathfile'], ['../bath/bath_1'])
        self.assertEqual(d['Defect']['axis_file'], 'Random3')   # the block's own spelling kept
        self.assertNotIn('defect', d)
        self.assertEqual(d['outfile'], 'name_1')
        self.assertIs(d['cluster_only'], True)
        self.assertEqual(d['method'], 'cCCE')

    def test_plain_copy(self):
        p = self.run_tool()
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(self.out) as f:
            self.assertEqual(json.load(f)['outfile'], 'run#1')

    def test_missing_block(self):
        p = self.run_tool('pcce.partition_size=2')
        self.assertEqual(p.returncode, 2)
        self.assertIn('pcce is not a block', p.stderr)
        self.assertFalse(os.path.exists(self.out))

    def test_out_is_in(self):
        p = self.run_tool('s_range=[1, 1]', out=self.src)
        self.assertEqual(p.returncode, 2)
        self.assertIn('OUT is IN', p.stderr)

    def test_not_key_value(self):
        self.assertEqual(self.run_tool('s_range').returncode, 2)

    def test_null_removes(self):
        """null leaves a key out (A94: configuration.s_range=null is the ensemble)."""
        p = self.run_tool('outfile=null', 'missing_key=null')
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(self.out) as f:
            d = json.load(f)
        self.assertNotIn('outfile', d)
        self.assertNotIn('missing_key', d)


if __name__ == '__main__':
    unittest.main()
