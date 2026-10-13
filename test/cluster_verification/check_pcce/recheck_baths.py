#!/usr/bin/env python3
"""Recheck immutable bath results and run only missing physical references.

Use --reference-results to recheck a downloaded recovery run without executing
the simulator. Original summaries, timings and outputs are never overwritten.
"""
import argparse
import collections
import copy
import hashlib
import json
from pathlib import Path
import re

from check_baths import compare, resolve, run_case, save, sha, values
from benchmark import parse_profiles


def execution(directory, dimension_bound=None):
    result = json.loads((directory/'execution.json').read_text())
    if result['returncode']:
        raise ValueError(f'Unsuccessful execution: {directory}')
    log = (directory/'run.log').read_text()
    profiles = parse_profiles(log)
    if not profiles:
        raise ValueError(f'Missing profiles: {directory}')
    config = json.loads((directory/'input.json').read_text())
    # the [profile] lines name the method as the code does: pcce, cce or gcce (A47: method cCCE/gCCE + cluster_mode)
    internal = 'pcce' if config.get('cluster_mode') == 'pcce' else 'gcce' if config['method'].lower() == 'gcce' else 'cce'
    if any(p['method'] != internal or
           p['backend'] != config['backend'] for p in profiles):
        raise ValueError(f'Unexpected method/backend: {directory}')
    if dimension_bound and max(int(p['dim']) for p in profiles) > dimension_bound:
        raise ValueError(f'Unexpected matrix dimension: {directory}')
    counts = collections.Counter()
    for profile in profiles:
        counts[profile['dim']] += int(profile['clusters'])
    clusters = sorted(line.strip() for line in log.splitlines()
                      if re.match(r'\s*Cluster\[\d+\]\[', line))
    if not clusters:
        raise ValueError(f'Missing cluster listing: {directory}')
    result.update(profiles=profiles, cluster_counts_by_dimension=dict(counts),
        cluster_listing_sha256=hashlib.sha256('\n'.join(clusters).encode()).hexdigest())
    compare(directory, directory, 1e-8, 2e-6)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original', type=Path, required=True)
    parser.add_argument('--fixtures', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cpu', type=Path)
    parser.add_argument('--reference-results', type=Path)
    args = parser.parse_args()
    args.output = args.output.resolve()
    args.fixtures = args.fixtures.resolve()
    original = json.loads((args.original/'summary.json').read_text())
    manifest = json.loads((args.fixtures/'manifest.json').read_text())
    if sha(args.fixtures/'manifest.json') != original['fixture_manifest_sha256']:
        raise ValueError('Changed fixture manifest')
    for name, expected in manifest['fixture_sha256'].items():
        if sha(args.fixtures/name) != expected:
            raise ValueError(f'Changed fixture: {name}')
    if args.cpu and sha(args.cpu) != original['binary_sha256']['cpu']:
        raise ValueError('CPU binary differs from original run')
    cases = {c['name']: c for c in manifest['cases']}
    if len(cases) != len(original['cases']) or set(cases) != {c['name'] for c in original['cases']}:
        raise ValueError('Incomplete or duplicate case coverage')
    args.output.mkdir(parents=True, exist_ok=False)
    summary = copy.deepcopy(original)
    summary.update(status='running', original_summary_sha256=sha(args.original/'summary.json'),
        original_status=original['status'], recovery='Reparse MPI profiles; preserve all original numerical runs.',
        additional_reference_runs=0)
    for record in summary['cases']:
        case = cases[record['name']]
        directory = args.original/case['name']
        record['original_status'] = record['status']
        if 'error' in record:
            record['original_error'] = record.pop('error')
        record.update(executions={}, comparisons={}, status='running')
        try:
            for backend in (['baseline'] if case['compare_baseline'] else [])+['cpu', 'cuda']:
                record['executions'][backend] = execution(directory/backend, case['maximum_dimension_bound'])
            for before, after in [('baseline', 'cpu'), ('cpu', 'cuda')]:
                if before not in record['executions']:
                    continue
                a, b = record['executions'][before], record['executions'][after]
                for key in ('cluster_listing_sha256', 'cluster_counts_by_dimension'):
                    if a[key] != b[key]:
                        raise ValueError(f'{before}/{after} {key} differs')
                record['comparisons'][f'{before}_{after}'] = compare(directory/before, directory/after, 1e-8, 2e-6)
            if case['reference']:
                reference = directory/'physical_reference'
                if not reference.exists():
                    reference = (args.reference_results or args.output)/case['name']/'physical_reference'
                    if not reference.exists():
                        if not args.cpu:
                            raise ValueError('Missing physical reference; supply --cpu to compute it')
                        run_case(dict(resolve(case['reference'], args.fixtures), backend='cpu'),
                                 reference, args.cpu.resolve(), original['ranks'], 240)
                        summary['additional_reference_runs'] += 1
                record['executions']['physical_reference'] = execution(reference)
                record['comparisons']['physical_reference'] = compare(reference, directory/'cpu', 1e-8, 2e-6)
            record['signal_change_from_initial'] = max(abs(row[1]-rows[0][1])
                for path in (directory/'cpu').glob('*_wiDiv') for rows in [values(path)] for row in rows)
            record['status'] = 'passed'
            print('PASS', case['name'], flush=True)
        except Exception as error:
            record.update(status='failed', error=f'{type(error).__name__}: {error}')
            print('FAIL', case['name'], record['error'], flush=True)
        save(args.output/'summary.json', summary)
    summary['status'] = 'passed' if all(c['status']=='passed' for c in summary['cases']) else 'failed'
    save(args.output/'summary.json', summary)
    if summary['status'] != 'passed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
