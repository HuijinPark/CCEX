"""Reject missing work and stale interpretations of bounded MPI diagnostics."""
import unittest
from partition_log import check_search, check_cache_hit


LOG = '''[partition-iteration] rank=0 trial=0
[partition-iteration] rank=0 trial=2
[partition-iteration] rank=0 trial=4
[partition] mode=deterministic total_trials=5 best_trial=3 best_rank=1
[partition-search] converged=2 cycles=2 limit=1 iterations=30 executed=20 skipped=10
'''


class PartitionLogTest(unittest.TestCase):
    def test_root_only(self):
        self.assertEqual(check_search(LOG, 5, 2)['total_trials'], 5)

    def test_bad_work_rejected(self):
        for before, after in [('total_trials=5', 'total_trials=4'),
                              ('cycles=2', 'cycles=1'),
                              ('executed=20', 'executed=19'),
                              ('best_rank=1', 'best_rank=0'),
                              ('rank=0 trial=2', 'rank=0 trial=0'),
                              ('rank=0 trial=4', 'rank=1 trial=4'),
                              ('[partition-iteration]', '[partition-trial]')]:
            with self.subTest(before=before), self.assertRaises(AssertionError):
                check_search(LOG.replace(before, after), 5, 2)

    def test_capped_details(self):
        capped = LOG.replace('[partition-iteration] rank=0 trial=4\n', '')
        capped += '[partition-detail] rank=0 suppressed=1 limit=2\n'
        self.assertEqual(len(check_search(capped, 5, 2, 2)['root_details']), 2)

    def test_cache_hit_has_no_search(self):
        hit = '[partition] mode=deterministic source=cache total_trials=0'
        check_cache_hit(hit)
        with self.assertRaises(AssertionError): check_cache_hit(hit + LOG)


if __name__ == '__main__':
    unittest.main()
