#!/usr/bin/env python3
"""Regressions for MPI profile forwarding and non-finite CCEX output."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'gpu_verification/check_gpu'))
from benchmark import parse_profiles, values


class ProfileParserTest(unittest.TestCase):
    def test_mpi_prefix_and_multiple_records(self):
        profile = ('[profile] rank=1 method=pcce backend=cpu dim=4 clusters=1 '
                   'build_s=0.000018355 setup_s=0.000041517 evolution_s=0.000038516')
        prefixed = '\t\t$items[  2][       3,5]'+profile
        records = parse_profiles(prefixed+'\n'+profile.replace('rank=1', 'rank=0'))
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0]['dim'], '4')
        self.assertEqual(sum(int(p['clusters']) for p in records), 2)
        self.assertEqual(records[1]['rank'], '0')

    def test_incomplete_record_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Incomplete'):
            parse_profiles('[profile] rank=0 method=pcce backend=cpu dim=4')


class CoherenceParserTest(unittest.TestCase):
    def test_nonfinite_output_including_fixed_width_padding(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'coherence_wiDiv'
            for token in ['+inf-inf         j', '-nan+inf         j',
                          '+inf-nanj', '-nan-nanj', '1e999+0j']:
                with self.subTest(token=token):
                    path.write_text('0\t1+0j\n0.0182547\t'+token+'\n')
                    with self.assertRaisesRegex(ValueError, 'Non-finite output: .*:2'):
                        values(path)

    def test_finite_complex_matrix_continuation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'coherence_wiDiv'
            path.write_text('0\t1+0j\t-0.5+0.2j\n  0.3-0.4j\n'
                            '0.1\t0+0j\t1e-30-2e-30j\t0.2+0j\n')
            self.assertEqual(values(path), [
                [0j, 1+0j, -.5+.2j, .3-.4j],
                [.1+0j, 0j, 1e-30-2e-30j, .2+0j]])


if __name__ == '__main__':
    unittest.main()
