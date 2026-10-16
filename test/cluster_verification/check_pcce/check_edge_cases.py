#!/usr/bin/env python3
"""pCCE empty-order regressions, with full-CCE/first-order ensemble references."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'gpu_verification/check_gpu'))
from benchmark import REPO, compare, put
from check_pcce import run_case


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpu', type=Path, required=True)
    parser.add_argument('--cuda', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--ranks', type=int, default=1)
    args = parser.parse_args()
    args.cpu = args.cpu.resolve()
    if args.cuda:
        args.cuda = args.cuda.resolve()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    original = json.loads((REPO/'test/gpu_verification/check_gpu/examples/pcce_ensemble.json').read_text())
    records = []
    for name, overrides in [('one-center', {'pcce.partition_size': 6}),
                            ('disconnected', {'rdip': 0.01}),
                            ('order4', {'order': 4})]:
        for state in (0, 1):
            config = put(original, dict(overrides, **{'pulse.count': 4, 'configuration.nstate': state}))
            directory = args.output/f'{name}_state{state}'
            cpu, _ = run_case(dict(config, backend='cpu'), directory/'cpu', args.cpu, args.ranks)
            assert cpu['returncode'] == 0 and cpu['profiles'], name
            compare(directory/'cpu', directory/'cpu', 1e-8, 2e-6)
            record = dict(case=directory.name, cpu=cpu)
            if state == 0:
                # Same physical full cluster, or independent partition groups.
                # Sampled states consume RNG during partitioning, so that
                # cross-method reference applies only to ensemble here.
                reference = dict(config, backend='cpu')
                if name == 'disconnected':
                    reference['order'] = 1
                else:
                    # the whole bath as one cCCE cluster: cluster_mode standard, without the pcce block
                    reference = put(reference, {'cluster_mode': 'standard', 'pcce': None, 'order': 6})
                result, _ = run_case(reference, directory/'reference', args.cpu, args.ranks)
                assert result['returncode'] == 0
                record['ensemble_reference'] = compare(directory/'reference', directory/'cpu', 1e-8, 2e-6)
            if args.cuda:
                gpu, _ = run_case(dict(config, backend='cuda'), directory/'cuda', args.cuda, args.ranks)
                assert gpu['returncode'] == 0 and gpu['profiles'], name
                assert {p['backend'] for p in gpu['profiles']} == {'cuda'}
                assert gpu['cluster_listing_sha256'] == cpu['cluster_listing_sha256']
                record['cuda'] = gpu
                record['cpu_cuda'] = compare(directory/'cpu', directory/'cuda', 1e-8, 2e-6)
            records.append(record)
            print('PASS', directory.name, flush=True)
    (args.output/'summary.json').write_text(json.dumps(dict(status='passed', ranks=args.ranks, cases=records), indent=2)+'\n')


if __name__ == '__main__':
    main()
