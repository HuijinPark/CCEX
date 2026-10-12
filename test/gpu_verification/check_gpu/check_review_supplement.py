#!/usr/bin/env python3
"""Generated-seed cache replay and actual large-vector batch admission."""
import argparse
import json
from pathlib import Path
import re
from benchmark import REPO, compare
from check_pcce import run_case


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cpu', type=Path, required=True)
    p.add_argument('--cuda', type=Path, required=True)
    p.add_argument('--before-cuda', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    config = json.loads((REPO/'test/gpu_verification/check_gpu/examples/pcce_ensemble.json').read_text())
    config.pop('seed')
    config.update(backend='cpu', pcce_cache_mode='readwrite',
                  pcce_cache_file=str(args.output/'partition.cache'))
    execution, log = run_case(config, args.output/'cache-generated', args.cpu.resolve(), ranks=2)
    assert execution['returncode'] == 0 and '[partition-cache] status=stored' in log, log[-2000:]
    seed = int(re.search(r'\[seed\].*shared with all ranks: (\d+);', log)[1])
    config.update(seed=seed, pcce_cache_mode='read')
    for ranks in (1, 4):
        directory = args.output/f'cache-replay-{ranks}'
        execution, log = run_case(config, directory, args.cpu.resolve(), ranks=ranks)
        assert execution['returncode'] == 0 and '[partition-cache] status=hit' in log, log[-2000:]
        compare(args.output/'cache-generated', directory, 1e-10, 1e-10)
    print('PASS generated-seed cache replay at 1/4 ranks', flush=True)
    records = []
    for method, fixture, order, dimension in (
        ('cce', 'p1_cce5_single.json', 4, 256),
        ('gcce', 'p1_gcce4_single.json', 3, 192),
    ):
        config = json.loads((REPO/'test/gpu_verification/check_gpu/examples'/fixture).read_text())
        config.update(order=order, nstep=64, npulse=4, gpu_batch_steps=64,
                      gpu_batch_clusters=8, gpu_cluster_batch_mb=64)
        entry = dict(method=method, dimension=dimension, order=order,
                     rbath=config['rbath'], rdip=config['rdip'], nstate=1, pulses=4,
                     steps=64, batch_budget_mib=64, executions={})
        for name, binary in (('cpu', args.cpu), ('before-cuda', args.before_cuda), ('cuda', args.cuda)):
            directory = args.output/method/name
            execution, log = run_case(dict(config, backend='cpu' if name=='cpu' else 'cuda'),
                                      directory, binary.resolve(), ranks=1, timeout=300)
            assert execution['returncode'] == 0, log[-2000:]
            error = compare(args.output/method/'cpu', directory, 1e-8, 2e-6)
            batches = [line for line in log.splitlines()
                       if '[gpu-batch]' in line and re.search(rf'\bdim={dimension}\b', line)]
            if name == 'before-cuda':
                assert not batches, batches
            if name == 'cuda':
                assert batches, f'{method} large-vector batches absent'
            entry['executions'][name] = dict(execution=execution, comparison=error, batches=batches)
            print('PASS vector admission', method, name, error, batches, flush=True)
        records.append(entry)
    (args.output/'summary.json').write_text(json.dumps(dict(status='passed', generated_seed=seed,
        cache_replay_ranks=[1, 4], cases=records), indent=2)+'\n')


if __name__ == '__main__':
    main()
