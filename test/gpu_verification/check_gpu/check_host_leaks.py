#!/usr/bin/env python3
"""Run normal-exit lifetime regressions with no CCEX leak suppressions.

Use --sanitize with a fully ASan/UBSan-instrumented CPU binary. With --cuda,
also compare the same inputs against CUDA single and cluster-batch execution.
Every output includes the exact input, log, and execution record.
"""
import argparse
import json
import os
from pathlib import Path

from benchmark import REPO, compare, get, put
from check_pcce import run_case


def cases():
    for method, fixture in [('cce', 'p1_cce5_single.json'), ('gcce', 'p1_gcce4_single.json')]:
        config = json.loads((REPO/'test/gpu_verification/check_gpu/examples'/fixture).read_text())
        config = put(config, {'order': 2, 'nstep': 9, 'pulse.count': 4, 'gpu_options.evolution_batch_steps': 4})
        for state in (0, 1, 3):
            yield f'p1-{method}-state{state}', put(config, {'configuration.nstate': state})
        if method == 'gcce':
            # One side's state given as a vector, the other from its m_s (Qubit.alpha_state/beta_state, A47 and A93;
            # an unpaired psia/psib was ignored before): its parsed temporary must still be released.
            for key in ('alpha_state', 'beta_state'):
                yield f'p1-gcce-unpaired-{key}', dict(config, Qubit=dict(config['Qubit'], **{key: {'vector': [1, 0, 0]}}))
    base = json.loads((REPO/'test/gpu_verification/check_gpu/examples/pcce_ensemble.json').read_text())
    for state in (0, 1):
        yield f'13c-pcce-state{state}', put(base, {'configuration.nstate': state, 'pulse.count': 4})
    # No Defect section, inline Qubit, and empty higher cluster orders.
    for method in ('cCCE', 'gCCE'):
        config = put(base, {'method': method, 'cluster_mode': 'standard', 'pcce': None, 'order': 4, 'rdip': 0.01})
        yield f'13c-{method.lower()}-empty-orders', config


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--cuda', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--ranks', type=int, default=1)
    p.add_argument('--sanitize', action='store_true')
    args = p.parse_args()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    if args.sanitize:
        os.environ['ASAN_OPTIONS'] = 'detect_leaks=1:halt_on_error=1'
        os.environ['UBSAN_OPTIONS'] = 'halt_on_error=1:print_stacktrace=1'
        os.environ.pop('LSAN_OPTIONS', None)
    summary = dict(status='running', ranks=args.ranks, sanitize=args.sanitize, cases=[])
    try:
        for name, config in cases():
            record = dict(case=name, executions={})
            variants = [('baseline', args.baseline, dict(backend='cpu')),
                        ('cpu', args.binary, dict(backend='cpu'))]
            if args.cuda:
                variants += [('cuda-single', args.cuda, dict(backend='cuda', gpu_options=dict(get(config, 'gpu_options', {}), evolution_batch_size=1))),
                             ('cuda-batch', args.cuda, dict(backend='cuda', gpu_options=dict(get(config, 'gpu_options', {}), evolution_batch_size=8)))]
            for variant, binary, overrides in variants:
                directory = args.output/name/variant
                execution, log = run_case(dict(config, **overrides), directory,
                                          binary.resolve(), ranks=args.ranks, timeout=180)
                if execution['returncode'] != 0:
                    raise RuntimeError(f'{name}/{variant} failed: {directory / "run.log"}')
                if any(marker in log for marker in ('ERROR: AddressSanitizer', 'ERROR: LeakSanitizer', 'runtime error:')):
                    raise RuntimeError(f'{name}/{variant} sanitizer diagnostic')
                atol, rtol = (1e-8, 2e-6) if args.sanitize or variant.startswith('cuda') else (1e-10, 1e-10)
                error = compare(args.output/name/'baseline', directory, atol, rtol)
                record['executions'][variant] = dict(**execution, comparison=error)
                print('PASS', name, variant, error, flush=True)
            summary['cases'].append(record)
        summary['status'] = 'passed'
    except BaseException as error:
        summary.update(status='failed', error=str(error))
        raise
    finally:
        (args.output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
