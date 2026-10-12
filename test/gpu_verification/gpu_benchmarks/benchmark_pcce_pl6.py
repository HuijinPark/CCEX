#!/usr/bin/env python3
"""Compare the supplied PL6 config-1 archive with same-server CPU/CUDA pCCE runs."""
import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import re

import numpy as np

from benchmark import REPO, values
from check_pcce import run_case


def metrics(reference, target, atol=1e-8, rtol=2e-6):
    if len(reference) != len(target):
        raise ValueError('Different output lengths')
    if not all(len(a) == len(b) and a[0] == b[0] for a, b in zip(reference, target)):
        raise ValueError('Different time grids or shapes')
    pairs = [(x, y) for a, b in zip(reference, target) for x, y in zip(a[1:], b[1:])]
    errors = [abs(x-y) for x, y in pairs]
    failures = [e > atol+rtol*abs(x) for e, (x, y) in zip(errors, pairs)]
    return dict(passed=not any(failures), values=len(pairs), failed_values=sum(failures),
                max_abs_error=max(errors), rms_error=float(np.sqrt(np.mean(np.square(errors)))))


def correction_function():
    # The repository's exact correction rule: coherence.py needs only numpy.
    source = REPO/'example/analysis/Analyzer/CoherenceAnalyzer.v1/coherence.py'
    spec = importlib.util.spec_from_file_location('coherence', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.error_correct, hashlib.sha256(source.read_bytes()).hexdigest()


def corrected(directory, correction):
    divided = values(directory/'coherence_wiDiv')
    undivided = values(directory/'coherence_noDiv')
    with contextlib.redirect_stdout(io.StringIO()):
        result = correction(np.array([row[1] for row in divided]),
                            np.array([row[1] for row in undivided]), 2.0)
    if not np.isfinite(result).all():
        raise ValueError('Non-finite corrected output')
    rows = [[row[0], complex(value)] for row, value in zip(divided, result)]
    (directory/'coherence_ErrCorr').write_text(''.join(
        f'{row[0].real:.10f}\t{row[1].real:+.10f}{row[1].imag:+.10f}j\n' for row in rows))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--ranks', type=int, default=32)
    parser.add_argument('--case', help='Optional input filename substring')
    parser.add_argument('--timeout', type=float, default=1800)
    parser.add_argument('--backends', nargs='+', choices=('cpu', 'cuda'), default=['cpu', 'cuda'])
    args = parser.parse_args()
    if args.backends not in (['cpu'], ['cpu', 'cuda']):
        parser.error('--backends must be cpu, or cpu cuda')
    args.package, args.output = args.package.resolve(), args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    correction, correction_sha = correction_function()
    summary = dict(status='running', spatial_config=1, ranks=args.ranks, repeats=1,
                   gpu_batch_steps=64, atol=1e-8, rtol=2e-6, cases=[], backends=args.backends,
                   correction_source_sha256=correction_sha,
                   note='Archived timings are historical log times on different hardware/rank counts; speedup uses new same-server runs only.')

    def save():
        tmp = args.output/'summary.json.tmp'
        tmp.write_text(json.dumps(summary, indent=2)+'\n')
        tmp.replace(args.output/'summary.json')

    save()
    try:
        for source in sorted((args.package/'inputs/ccein').glob('*.json')):
            if args.case and args.case not in source.name:
                continue
            config = json.loads(source.read_text())
            grid = 'coarse_grid_nstep500' if 'coarse' in source.name else 'fine_grid'
            name = re.sub(r'^ccein_(coarse|fine)_', '', source.stem)
            historical = args.package/'results'/grid/name
            old_log = (historical/'process_1').read_text()
            seed = int(re.search(r'^\s*seed\s*:\s*(\d+)', old_log, re.M)[1])
            config['gyrofile'] = str((source.parent/config['gyrofile']).resolve())
            config['bathfile'] = [str((source.parent/p).resolve()) for p in config['bathfile']]
            config.update(seed=seed, profile=True, savemode='normal', gpu_options=dict(config.get('gpu_options', {}), evolution_batch_steps=64))
            directory = args.output/source.stem
            old = directory/'archived'
            old.mkdir(parents=True)
            for suffix in ('noDiv', 'wiDiv', 'ErrCorr'):
                matches = list(historical.glob(f'*_config1_{suffix}'))
                if len(matches) != 1:
                    raise ValueError(f'Expected one archived config-1 {suffix}: {historical}')
                (old/f'coherence_{suffix}').write_bytes(matches[0].read_bytes())
            record = dict(case=source.stem, seed=seed, sK=config['sK'], npulse=config['npulse'],
                          nstep=config['nstep'], rdip=config['rdip'], backends={},
                          archived_wall_s=float(re.search(r'JOB DONE : Wall time = ([\d.]+)', old_log)[1]))
            summary['cases'].append(record)
            summary['current_case'] = source.stem
            save()
            archive_correction = corrected(old, correction)
            # Preserve/check the original archive correction before regenerating it.
            original_corrected = values(next(historical.glob('*_config1_ErrCorr')))
            record['archive_correction_recheck'] = metrics(original_corrected, archive_correction)
            for backend in args.backends:
                summary['current_backend'] = backend
                save()
                result, log = run_case(dict(config, backend=backend), directory/backend,
                                       REPO/f'bin/main.{backend}.out', args.ranks, args.timeout)
                record['backends'][backend] = result
                save()
                if result['returncode'] or not result['profiles']:
                    raise RuntimeError(f'{source.stem} {backend} failed; inspect run.log')
                expected = 'cuda' if backend == 'cuda' else 'cpu'
                if any(p['method'] != 'pcce' or p['backend'] != expected for p in result['profiles']):
                    raise ValueError('Unexpected pCCE backend profile')
                result['against_archive'] = {suffix: metrics(values(old/f'coherence_{suffix}'),
                    values(directory/backend/f'coherence_{suffix}')) for suffix in ('noDiv', 'wiDiv')}
                new_corrected = corrected(directory/backend, correction)
                result['corrected_against_archive'] = metrics(original_corrected, new_corrected)
                if backend == 'cuda':
                    cpu = record['backends']['cpu']
                    record['same_cluster_partition'] = cpu['cluster_listing_sha256'] == result['cluster_listing_sha256']
                    record['cpu_cuda'] = {suffix: metrics(values(directory/'cpu'/f'coherence_{suffix}'),
                        values(directory/'cuda'/f'coherence_{suffix}')) for suffix in ('noDiv', 'wiDiv', 'ErrCorr')}
                    record['speedup'] = cpu['wall_s']/result['wall_s']
                save()
                print(source.stem, backend, 'wall_s=', round(result['wall_s'], 3),
                      'archive=', {k:v['passed'] for k,v in result['against_archive'].items()}, flush=True)
            if 'cuda' in args.backends:
                print('CPU/CUDA', record['cpu_cuda'], 'speedup=', record['speedup'], flush=True)
        if not summary['cases']:
            raise ValueError('No matching inputs')
        if 'cuda' in args.backends:
            paired = all(c['same_cluster_partition'] and all(v['passed'] for v in c['cpu_cuda'].values()) for c in summary['cases'])
            summary['status'] = 'passed' if paired else 'comparison_failed'
        else:
            summary['status'] = 'cpu_reference_completed'
        summary['all_archive_comparisons_passed'] = all(c['archive_correction_recheck']['passed'] and
            all(v['passed'] for b in c['backends'].values() for v in b['against_archive'].values()) and
            all(b['corrected_against_archive']['passed'] for b in c['backends'].values()) for c in summary['cases'])
        save()
    except BaseException as error:
        summary.update(status='failed', error=f'{type(error).__name__}: {error}')
        save()
        raise
    if summary['status'] == 'comparison_failed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
