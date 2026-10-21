#!/usr/bin/env python3
"""Repeat selected bath cases with alternating CPU/CUDA execution order."""
import argparse
import json
from pathlib import Path
import statistics

from check_baths import compare, resolve, run_case, save, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixtures', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cpu', type=Path, required=True)
    parser.add_argument('--cuda', type=Path, required=True)
    parser.add_argument('--case', action='append', required=True, help='Exact case name; repeat for more cases')
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 2 or len(set(args.case)) != len(args.case):
        parser.error('At least two repeats and distinct case names are required')
    args.output = args.output.resolve()
    args.fixtures = args.fixtures.resolve()
    manifest = json.loads((args.fixtures/'manifest.json').read_text())
    original = json.loads((args.reference/'summary.json').read_text())
    if sha(args.fixtures/'manifest.json') != original['fixture_manifest_sha256']:
        raise ValueError('Changed fixture manifest')
    for name, expected in manifest['fixture_sha256'].items():
        if sha(args.fixtures/name) != expected:
            raise ValueError(f'Changed fixture: {name}')
    binaries = dict(cpu=args.cpu.resolve(), cuda=args.cuda.resolve())
    for name, path in binaries.items():
        if sha(path) != original['binary_sha256'][name]:
            raise ValueError(f'Changed {name} binary')
    cases = {c['name']: c for c in manifest['cases']}
    if not all(name in cases for name in args.case):
        parser.error('Unknown case')
    args.output.mkdir(parents=True, exist_ok=False)
    summary = dict(status='running', repeats=args.repeats, ranks=original['ranks'],
        fixture_manifest_sha256=sha(args.fixtures/'manifest.json'),
        binary_sha256={name: sha(path) for name, path in binaries.items()},
        timing='Fresh process wall time including MPI startup, clustering, CUDA initialization and output; no warmup; alternating backend order.',
        cases=[])
    for index, name in enumerate(args.case):
        case = cases[name]
        record = dict(name=name, trials=[])
        summary['cases'].append(record)
        for repeat in range(args.repeats):
            trial = dict(index=repeat, executions={}, comparisons={})
            record['trials'].append(trial)
            order = ['cpu', 'cuda'] if (repeat+index)%2 == 0 else ['cuda', 'cpu']
            trial['order'] = order
            root = args.output/name/f'repeat{repeat}'
            for backend in order:
                result, _ = run_case(dict(resolve(case['config'], args.fixtures), backend=backend),
                    root/backend, binaries[backend], original['ranks'], 240)
                if result['returncode'] or not result['profiles'] or any(
                        p['backend'] != backend or p['method'] != 'pcce' for p in result['profiles']):
                    raise ValueError(f'Failed {name} {backend} execution')
                if max(int(p['dim']) for p in result['profiles']) > case['maximum_dimension_bound']:
                    raise ValueError('Unexpected matrix dimension')
                trial['executions'][backend] = result
                trial['comparisons'][backend] = compare(args.reference/name/'cpu', root/backend, 1e-8, 2e-6)
                save(args.output/'summary.json', summary)
            for key in ('cluster_counts_by_dimension', 'cluster_listing_sha256'):
                if trial['executions']['cpu'][key] != trial['executions']['cuda'][key]:
                    raise ValueError(f'Changed partition: {name}')
            print('PASS', name, repeat, *(f'{b}={trial["executions"][b]["wall_s"]:.3f}s' for b in ('cpu', 'cuda')), flush=True)
        record['median_s'] = {b: statistics.median(t['executions'][b]['wall_s'] for t in record['trials']) for b in ('cpu', 'cuda')}
        record['cpu_over_cuda'] = record['median_s']['cpu']/record['median_s']['cuda']
        save(args.output/'summary.json', summary)
    summary['status'] = 'passed'
    save(args.output/'summary.json', summary)


if __name__ == '__main__':
    main()
