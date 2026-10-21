"""Check bounded root-only details and MPI-reduced pCCE search counters."""
import re


def check_search(log, trials, ranks, detail_limit=32):
    summaries = re.findall(
        r'\[partition\] mode=deterministic total_trials=(\d+) best_trial=(\d+) best_rank=(\d+)', log)
    assert len(summaries) == 1, summaries
    total, best, owner = map(int, summaries[0])
    assert total == trials and 0 <= best < trials and owner == best % ranks, summaries
    counters = re.findall(
        r'\[partition-search\] converged=(\d+) cycles=(\d+) limit=(\d+) '
        r'iterations=(\d+) executed=(\d+) skipped=(\d+)', log)
    assert len(counters) == 1, counters
    stats = dict(zip(('converged', 'cycles', 'limit', 'iterations', 'executed', 'skipped'),
                     map(int, counters[0])))
    assert stats['converged'] + stats['cycles'] + stats['limit'] == trials, stats
    assert stats['executed'] + stats['skipped'] == stats['iterations'], stats
    details = [(int(rank), int(trial)) for rank, trial in re.findall(
        r'\[partition-iteration\] rank=(\d+) trial=(\d+)', log)]
    # Other ranks deliberately do not print per-trial details. Validate their
    # reduced totals above; do not mistake the bounded root log for all trials.
    expected = [(0, trial) for trial in range(0, trials, ranks)][:detail_limit]
    assert details == expected, (details, expected)
    suppressed = re.findall(r'\[partition-detail\] rank=0 suppressed=(\d+) limit=(\d+)', log)
    remaining = len(range(0, trials, ranks)) - len(expected)
    assert suppressed == ([(str(remaining), str(detail_limit))] if remaining else []), suppressed
    return dict(total_trials=total, best_trial=best, best_rank=owner, root_details=details, search=stats)


def check_cache_hit(log):
    assert '[partition-iteration]' not in log and '[partition-search]' not in log
    assert 'source=cache total_trials=0' in log
