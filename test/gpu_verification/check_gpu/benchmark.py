#!/usr/bin/env python3
"""Run immutable regression inputs in an isolated output directory, using stdlib only.

The reference uses backend=cpu; repeated runs alternate execution order.
Numeric comparisons include both raw coherence products, with finite-value and
matching time-grid checks.
"""
import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess
import statistics
import time

REPO = Path(__file__).resolve().parents[3]
INPUTS = REPO / 'test/code_verification/check_input_forms/inputs'


def parse_profiles(log):
    # MPI's output forwarding can concatenate another rank's unterminated
    # output with a complete profile record. A profile need not start a line.
    number = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
    pattern = (r'\[profile\]\s+(rank=\d+\s+method=\w+\s+backend=\w+\s+'
               r'dim=\d+\s+clusters=\d+\s+build_s='+number+
               r'\s+setup_s='+number+r'\s+evolution_s='+number+')')
    records = [dict(part.split('=', 1) for part in match.split())
               for match in re.findall(pattern, log)]
    if len(records) != log.count('[profile]'):
        raise ValueError('Incomplete or interleaved profile record; inspect run.log')
    return records


def values(path):
    rows = []
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        parts = line.split()
        if not parts:
            continue
        # CCEX's fixed-width %lf output can put padding between an infinite
        # imaginary component and 'j'. Report the numerical failure before
        # splitting that token, rather than hiding it as malformed complex text.
        if re.search(r'nan|inf', line, re.IGNORECASE):
            raise ValueError(f'Non-finite output: {path}:{line_number}')
        row = [complex(x) for x in parts]
        if any(not math.isfinite(v.real) or not math.isfinite(v.imag) for v in row):
            raise ValueError(f'Non-finite output: {path}:{line_number}')
        if 'j' in parts[0]:
            if not rows:
                raise ValueError(f'Missing time before matrix continuation: {path}')
            rows[-1].extend(row)
        else:
            rows.append(row)
    if not rows:
        raise ValueError(f'Empty output: {path}')
    return rows


def compare(reference, target, atol, rtol):
    ref_files = sorted(reference.glob('*_noDiv')) + sorted(reference.glob('*_wiDiv'))
    new_files = sorted(target.glob('*_noDiv')) + sorted(target.glob('*_wiDiv'))
    if not ref_files or [p.name for p in ref_files] != [p.name for p in new_files]:
        raise ValueError('Missing or mismatched coherence files')
    maximum, squared, count = 0.0, 0.0, 0
    for ref, new in zip(ref_files, new_files):
        a, b = values(ref), values(new)
        if len(a) != len(b):
            raise ValueError('Different row counts')
        for ra, rb in zip(a, b):
            if len(ra) != len(rb) or ra[0] != rb[0]:
                raise ValueError('Different time grids/output shapes')
            for x, y in zip(ra[1:], rb[1:]):
                error = abs(x-y)
                if error > atol + rtol*abs(x):
                    raise ValueError(f'{new.name}: |CPU-target|={error:.3g} exceeds tolerance')
                maximum = max(maximum, error)
                squared += error*error
                count += 1
    return {'max_abs_error': maximum, 'rms_error': math.sqrt(squared/count), 'values': count}


def run(binary, source, directory, backend, args):
    directory.mkdir(parents=True, exist_ok=False)
    config = json.loads(source.read_text())
    for key in ['qubitfile', 'gyrofile', 'avaaxfile', 'statefile', 'exstatefile']:
        if config.get(key):
            config[key] = str((source.parent.parent / config[key]).resolve())
    config['bathfile'] = [str((source.parent.parent / p).resolve()) for p in config['bathfile']]
    if args.nearest_sites:
        if len(config['bathfile']) != 1:
            raise ValueError('--nearest-sites requires one bath file')
        if 'Qubit' in config:
            origin = config['Qubit']['qubit'][0]['xyz']
        else:
            origin = [float(v) for v in Path(config['qubitfile']).read_text().split()[:3]]
            if len(origin)!=3:
                raise ValueError('Cannot read qubit origin')
        rows = Path(config['bathfile'][0]).read_text().splitlines()[1:]
        ranked = sorted(enumerate(rows),key=lambda item: sum(
            (float(v)-c)**2 for v,c in zip(item[1].split()[:3],origin)))
        selected = sorted(ranked[:args.nearest_sites])
        if len(selected) != args.nearest_sites:
            raise ValueError('Not enough bath sites')
        bath = directory/'bath'
        bath.write_text(f'{len(selected)+1}.0 0 0 0\n'+'\n'.join(row for _,row in selected)+'\n')
        config['bathfile'] = [str(bath)]
        radius = max(math.dist([float(v) for v in row.split()[:3]],origin) for _,row in selected)+1
        config.update(order=args.nearest_sites,rbath=radius,rdip=2*radius)
    config.update(backend=backend, profile=True, savemode=args.savemode, outfile=str(directory/'coherence'))
    if args.pulses is not None:
        config['npulse'] = args.pulses
    if args.steps is not None:
        config['nstep'] = args.steps
    if args.batch is not None:
        config['gpu_batch_steps'] = args.batch
    if args.min_dim is not None:
        config['gpu_min_dim'] = args.min_dim
    if args.evolution is not None:
        config['evolution'] = args.evolution
    path = directory/'input.json'
    path.write_text(json.dumps(config, indent=2)+'\n')
    command = [str(binary), '-f', str(path)]
    if args.mpirun:
        command = [args.mpirun, '-n', str(args.ranks)] + command
    elif args.ranks != 1:
        raise ValueError('--ranks requires --mpirun')
    env = dict(os.environ)
    env.setdefault('I_MPI_FABRICS', 'shm')
    start = time.perf_counter()
    with (directory/'run.log').open('w') as log:
        subprocess.run(command, cwd=REPO, env=env, stdout=log, stderr=subprocess.STDOUT,
                       timeout=args.timeout, check=True)
    elapsed = time.perf_counter()-start
    output = (directory/'run.log').read_text()
    profiles = parse_profiles(output)
    if backend == 'cuda' and not any(p['backend'] == 'cuda' for p in profiles):
        raise ValueError('Requested CUDA but no GPU computation was profiled')
    return {'wall_s': elapsed, 'profiles': profiles}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--backend', choices=['cpu', 'cuda', 'auto'], default='cuda')
    parser.add_argument('--method', choices=['cce', 'gcce', 'all'], default='all')
    parser.add_argument('--case', help='Only case names containing this text')
    parser.add_argument('--nearest-sites', type=int, help='Controlled size experiment: keep N nearest bath sites and use order N')
    parser.add_argument('--output', type=Path, required=True, help='New directory; existing results are never overwritten')
    parser.add_argument('--steps', type=int)
    parser.add_argument('--savemode', choices=['normal','allfull'], default='normal')
    parser.add_argument('--pulses', type=int)
    parser.add_argument('--batch', type=int)
    parser.add_argument('--min-dim', type=int)
    parser.add_argument('--evolution', choices=['vector', 'matrix'])
    parser.add_argument('--nstate', choices=['0', '1', 'both'], default='both')
    parser.add_argument('--mpirun')
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--timeout', type=float, default=300)
    parser.add_argument('--repeats', type=int, default=1, help='Interleave CPU/target runs; report median wall times')
    parser.add_argument('--atol', type=float, default=1e-8)
    # Long-time NV/13C cluster division amplifies eigensolver rounding. The
    # independent CPU eigen-vs-expm comparison shows the same ppm-scale effect.
    # Small-system CUDA tests separately enforce 3e-11 against matrix exp().
    parser.add_argument('--rtol', type=float, default=2e-6)
    args = parser.parse_args()
    if args.repeats<1:
        parser.error('--repeats must be positive')
    if args.nearest_sites is not None and args.nearest_sites<1:
        parser.error('--nearest-sites must be positive')
    if args.backend=='cuda' and args.savemode=='allfull' and args.method!='gcce':
        parser.error('cCCE does not export full density matrices; use --method gcce for --savemode allfull')
    args.binary, args.reference, args.output = (p.resolve() for p in (args.binary,args.reference,args.output))
    args.output.mkdir(parents=True, exist_ok=False)
    summary = []
    for case in sorted(INPUTS.glob('Diamond_*')):
        if args.case and args.case not in case.name:
            continue
        method = case.name.rsplit('_',1)[1]
        if args.method != 'all' and args.method != method:
            continue
        for state in [0,1]:
            if args.nstate != 'both' and int(args.nstate) != state:
                continue
            source = case/'inputfiles'/f'ccein_nstate{state}.json'
            root = args.output/f'{case.name}_nstate{state}'
            baselines, targets, errors = [], [], []
            for repeat in range(args.repeats):
                directory=root if args.repeats==1 else root/f'repeat{repeat+1}'
                # Alternate execution order to limit warm-cache/order bias.
                if repeat%2:
                    target=run(args.binary,source,directory/'target',args.backend,args)
                    baseline=run(args.reference,source,directory/'reference','cpu',args)
                else:
                    baseline=run(args.reference,source,directory/'reference','cpu',args)
                    target=run(args.binary,source,directory/'target',args.backend,args)
                baselines.append(baseline); targets.append(target)
                errors.append(compare(directory/'reference',directory/'target',args.atol,args.rtol))
            baseline=dict(baselines[-1],wall_s=statistics.median(r['wall_s'] for r in baselines))
            target=dict(targets[-1],wall_s=statistics.median(r['wall_s'] for r in targets))
            error=dict(max_abs_error=max(e['max_abs_error'] for e in errors),
                       rms_error=max(e['rms_error'] for e in errors),values=errors[0]['values'])
            record = dict(case=root.name, reference=baseline, target=target, **error)
            record.update(repeats=args.repeats,reference_runs=baselines,target_runs=targets)
            record['speedup'] = baseline['wall_s']/target['wall_s']
            summary.append(record)
            (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
            print(f"{root.name}: max_error={error['max_abs_error']:.3g} "
                  f"CPU={baseline['wall_s']:.3f}s target={target['wall_s']:.3f}s "
                  f"speedup={record['speedup']:.3f}",flush=True)
    if not summary:
        raise ValueError('No test cases matched the requested method/case/state')


if __name__ == '__main__':
    main()
