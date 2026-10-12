#!/usr/bin/env python3
"""Small pCCE integration checks: partitions, pulse/state modes and CUDA selection."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

from benchmark import REPO, compare, parse_profiles


# mpirun as test/README.md has it: CCEX_MPI when set, else the Intel MPI of the loaded module, else PATH. With
# ~/anaconda3/bin first on PATH (as the Python checks need on WSL), a bare mpirun - and the mpiexec.hydra and proxy
# Intel's mpirun looks up on PATH - were anaconda's MPICH, which does not launch Intel MPI ranks as one job: each ran
# alone into the same run.log, and check_partition_cache lost refusal messages in about one run in two (2026-10-11).
# So the launcher's own directory goes first on the child's PATH too.
def mpirun():
    if os.environ.get('CCEX_MPI'):
        return os.environ['CCEX_MPI']
    intel = os.path.join(os.environ.get('I_MPI_ROOT', ''), 'bin', 'mpirun')
    return intel if os.environ.get('I_MPI_ROOT') and os.access(intel, os.X_OK) else 'mpirun'


def run_case(config, directory, binary, ranks=1, timeout=120, hidden=False, extra=()):
    directory.mkdir(parents=True, exist_ok=False)
    config = dict(config, outfile=str(directory / 'coherence'))
    path = directory / 'input.json'
    path.write_text(json.dumps(config, indent=2) + '\n')
    launcher = mpirun()
    command = [launcher, '-n', str(ranks),
               str(binary), '-f', str(path), '-v', *extra]
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    env.setdefault('I_MPI_FABRICS', 'shm')
    if os.path.dirname(launcher):
        env['PATH'] = os.path.dirname(launcher) + os.pathsep + env.get('PATH', '')
    if hidden:
        env['CUDA_VISIBLE_DEVICES'] = ''
    start = time.perf_counter()
    with (directory / 'run.log').open('w') as output:
        result = subprocess.run(command, cwd=REPO, env=env, stdout=output,
                                stderr=subprocess.STDOUT, timeout=timeout)
    wall = time.perf_counter() - start
    log = (directory / 'run.log').read_text()
    profiles = parse_profiles(log)
    clusters = sorted(line.strip() for line in log.splitlines()
                      if re.match(r'\s*Cluster\[\d+\]\[', line))
    counts = {}
    for profile in profiles:
        dim = profile['dim']
        counts[dim] = counts.get(dim, 0) + int(profile['clusters'])
    record = dict(wall_s=wall, returncode=result.returncode, profiles=profiles,
                  cluster_counts_by_dimension=counts,
                  cluster_listing_sha256=hashlib.sha256('\n'.join(clusters).encode()).hexdigest())
    (directory / 'execution.json').write_text(json.dumps(record, indent=2) + '\n')
    return record, log


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpu', type=Path, required=True)
    parser.add_argument('--cuda', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--ranks', type=int, default=1)
    args = parser.parse_args()
    args.cpu, args.cuda, args.output = (p.resolve() for p in (args.cpu, args.cuda, args.output))
    args.output.mkdir(parents=True, exist_ok=False)
    original = json.loads((REPO / 'test/gpu_verification/check_gpu/examples/pcce_ensemble.json').read_text())
    records = []
    for size in (1, 2, 3):
        for state in (0, 1):
            for pulses in (1, 4):
                name = f'sK{size}_state{state}_pulse{pulses}'
                config = dict(original, sK=size, nstate=state, npulse=pulses)
                cpu, _ = run_case(dict(config, backend='cpu'), args.output/name/'cpu', args.cpu, args.ranks)
                assert cpu['returncode'] == 0, name
                compare(args.output/name/'cpu', args.output/name/'cpu', 1e-8, 2e-6)
                for backend in ('cuda', 'auto'):
                    gpu, _ = run_case(dict(config, backend=backend), args.output/name/backend, args.cuda, args.ranks)
                    assert gpu['returncode'] == 0 and gpu['profiles'], (name, backend)
                    assert all(p['method'] == 'pcce' for p in gpu['profiles'])
                    assert gpu['cluster_counts_by_dimension'] == cpu['cluster_counts_by_dimension']
                    assert gpu['cluster_listing_sha256'] == cpu['cluster_listing_sha256']
                    expected = {'cuda'} if backend == 'cuda' else ({'cpu'} if size == 1 else {'cpu', 'cuda'} if size == 2 else {'cuda'})
                    assert {p['backend'] for p in gpu['profiles']} == expected, (name, backend, gpu['profiles'])
                    error = compare(args.output/name/'cpu', args.output/name/backend, 1e-8, 2e-6)
                    records.append(dict(case=name, backend=backend, cpu_s=cpu['wall_s'], cuda_s=gpu['wall_s'], **error))
                    print('PASS', name, backend, error, flush=True)

    for name, overrides, expected in [
        ('missing-sK', {'sK': None}, 'pCCE requires sK'),
        ('zero-sK', {'sK': 0}, 'positive 32-bit integer'),
        ('fractional-sK', {'sK': 1.5}, 'positive 32-bit integer'),
        ('string-sK', {'sK': '2'}, 'finite JSON number'),
        ('zero-trials', {'max_trial': 0}, 'positive 32-bit integer'),
        ('zero-iterations', {'max_iter': 0}, 'positive 32-bit integer'),
        ('oversized-sK', {'sK': 7}, 'bath spin count'),
    ]:
        config = dict(original, backend='cpu', **overrides)
        if config.get('sK') is None:
            config.pop('sK')
        result, log = run_case(config, args.output/name, args.cpu)
        assert result['returncode'] != 0 and expected in log, (name, log[-1000:])
        print('PASS', name, flush=True)

    reference = args.output/'sK2_state0_pulse1/cpu'
    for name, binary, backend, hidden, expected in [
        ('cpu-auto', args.cpu, 'auto', False, None),
        ('cpu-cuda-rejected', args.cpu, 'cuda', False, 'no CUDA support'),
        ('hidden-auto', args.cuda, 'auto', True, None),
        ('hidden-cuda-rejected', args.cuda, 'cuda', True, 'CUDA initialization failed'),
    ]:
        result, log = run_case(dict(original, backend=backend), args.output/name, binary, args.ranks, hidden=hidden)
        if expected:
            assert result['returncode'] != 0 and expected in log, name
        else:
            assert result['returncode'] == 0 and result['profiles'], name
            assert {p['backend'] for p in result['profiles']} == {'cpu'}, name
            compare(reference, args.output/name, 1e-8, 2e-6)
        print('PASS', name, flush=True)
    (args.output/'summary.json').write_text(json.dumps(dict(status='passed', ranks=args.ranks,
        comparisons=records, input_checks=7, fallback_checks=4), indent=2) + '\n')


if __name__ == '__main__':
    main()
