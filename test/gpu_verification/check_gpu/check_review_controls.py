#!/usr/bin/env python3
"""Replay a generated MPI pCCE seed and verify bounded rank-zero details."""
import argparse
import json
from pathlib import Path
import re
from benchmark import REPO, compare
from check_pcce import run_case
from check_seed_states import sampled, random_reference, check_stream


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpu', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.cpu, args.output = args.cpu.resolve(), args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    config = json.loads((REPO/'test/gpu_verification/check_gpu/examples/pcce_ensemble.json').read_text())
    config.update(backend='cpu', nstate=2, max_trial=12, iter_detail=True, pcce_detail_limit=2)
    config.pop('seed')
    first, log = run_case(config, args.output/'generated', args.cpu, ranks=2)
    assert first['returncode'] == 0, log[-2000:]
    seeds = re.findall(r'\[seed\] generated on rank 0 and shared with all ranks: (\d+);', log)
    assert len(seeds) == 1, seeds
    details = re.findall(r'\[partition-iteration\] rank=(\d+)', log)
    assert details == ['0', '0'], details
    config['seed'] = int(seeds[0])
    states = sampled(log)
    check_stream(states, random_reference(config['seed'], 6, 2))
    for ranks in (1, 2, 4):
        directory = args.output/f'replay-{ranks}'
        replay, output = run_case(dict(config, pcce_detail_limit=0), directory, args.cpu, ranks=ranks)
        assert replay['returncode'] == 0, output[-2000:]
        assert '[partition-iteration]' not in output
        assert replay['cluster_listing_sha256'] == first['cluster_listing_sha256']
        assert sampled(output) == states
        compare(args.output/'generated', directory, 1e-10, 1e-10)
    (args.output/'summary.json').write_text(json.dumps(dict(status='passed', generated_seed=int(seeds[0]),
        source_ranks=2, replay_ranks=[1, 2, 4], detail_lines=2), indent=2)+'\n')
    print('PASS generated MPI seed replay at 1/2/4 ranks and bounded rank-zero logging')


if __name__ == '__main__':
    main()
