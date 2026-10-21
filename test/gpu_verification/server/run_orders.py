#!/usr/bin/env python3
"""Compare orders 3--5 on a whole-server Slurm allocation, one order per job."""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'check_gpu'))  # benchmark.py
from benchmark import INPUTS, REPO, a94, compare, absolute_paths, get, put


def prepare_config(source, directory, order, backend, batch):
    config = json.loads(source.read_text())
    config = absolute_paths(config, source.parent.parent)
    config['order'] = order
    return put(config, {'backend': backend, 'profile': True, 'savemode': 'coherence', 'outfile': str(directory / 'coherence'),
                        'gpu_options': dict(get(config, 'gpu_options', {}), evolution_batch_steps=batch, min_dim=256)})


def sample_resources():
    memory = {}
    for line in Path('/proc/meminfo').read_text().splitlines():
        parts = line.split()
        if parts[0] in ('MemTotal:', 'MemAvailable:'):
            memory[parts[0][:-1]] = int(parts[1]) * 1024
    result = subprocess.run(
        ['nvidia-smi', '--query-gpu=memory.used,memory.total', '--format=csv,noheader,nounits'],
        capture_output=True, text=True, check=True, timeout=10)
    devices = [[int(v.strip()) for v in line.split(',')] for line in result.stdout.splitlines() if line.strip()]
    if len(devices) != 1:
        raise RuntimeError('This whole-server benchmark expects exactly one physical GPU')
    used, total = devices[0]
    return dict(host_available_bytes=memory['MemAvailable'], host_total_bytes=memory['MemTotal'],
                gpu_used_mib=used, gpu_total_mib=total)


def stop_process(process):
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        process.wait(timeout=10)
        return
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait(timeout=10)


def monitored_run(command, directory, timeout, sampler=sample_resources,
                  min_available_bytes=8 * 1024**3, max_gpu_fraction=0.90, interval=2):
    """Monitor whole-node memory; stop this case and let Slurm clean up on failure."""
    start = time.perf_counter()
    peak = dict(host_used_bytes=0, gpu_used_mib=0)
    with (directory / 'run.log').open('w') as log, (directory / 'resources.jsonl').open('w') as resources:
        process = subprocess.Popen(command, cwd=REPO, stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        try:
            while True:
                sample = sampler()
                sample['elapsed_s'] = time.perf_counter() - start
                resources.write(json.dumps(sample) + '\n')
                resources.flush()
                peak['host_used_bytes'] = max(peak['host_used_bytes'], sample['host_total_bytes'] - sample['host_available_bytes'])
                peak['gpu_used_mib'] = max(peak['gpu_used_mib'], sample['gpu_used_mib'])
                if sample['host_available_bytes'] < min_available_bytes:
                    raise MemoryError('Whole-node available RAM fell below the 8 GiB reserve')
                if sample['gpu_used_mib'] > max_gpu_fraction * sample['gpu_total_mib']:
                    raise MemoryError('GPU memory exceeded the 90 percent guard')
                remaining = timeout - (time.perf_counter() - start)
                if remaining <= 0:
                    raise TimeoutError('Per-case wall-time limit exceeded')
                try:
                    code = process.wait(timeout=min(interval, remaining))
                    if code:
                        raise RuntimeError(f'Calculation exited with status {code}; see {directory / "run.log"}')
                    break
                except subprocess.TimeoutExpired:
                    pass
        except BaseException:
            stop_process(process)
            raise
    return dict(wall_s=time.perf_counter() - start, sampled_peak=peak)


def run_case(source, directory, order, backend, args):
    directory.mkdir(parents=True, exist_ok=False)
    config = prepare_config(source, directory, order, backend, args.batch)
    path = directory / 'input.json'
    path.write_text(json.dumps(a94(config), indent=2) + '\n')   # the A94 forms
    binary = REPO / ('bin/main.out' if backend == 'cpu' else 'bin/main_gpu.out')
    command = ['bash', str(REPO / 'test/gpu_verification/server/mpi_launch.sh'), '-n', str(args.ranks),
               str(binary), '-f', str(path)]
    result = monitored_run(command, directory, args.timeout)
    result['profiles'] = [dict(item.split('=', 1) for item in line.split()[1:])
                          for line in (directory / 'run.log').read_text().splitlines()
                          if line.startswith('[profile] ')]
    if not result['profiles']:
        raise ValueError('Calculation returned no cluster profiles')
    if backend == 'cuda' and any(p['backend'] != 'cuda' for p in result['profiles']):
        raise ValueError('Forced CUDA did not use CUDA for every profiled cluster')
    result['cluster_counts_by_dimension'] = {}
    for p in result['profiles']:
        key = p['dim']
        result['cluster_counts_by_dimension'][key] = result['cluster_counts_by_dimension'].get(key, 0) + int(p['clusters'])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--order', type=int, choices=(3, 4, 5), required=True)
    parser.add_argument('--ranks', type=int, default=32)
    parser.add_argument('--batch', type=int, default=4)
    parser.add_argument('--timeout', type=float, default=72 * 3600)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.ranks < 1 or args.batch < 1 or args.timeout <= 0:
        parser.error('ranks, batch and timeout must be positive')
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Run this driver inside the Slurm allocation')

    def terminate(signum, frame):
        raise SystemExit(f'Received signal {signum}')

    signal.signal(signal.SIGTERM, terminate)
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    summary = dict(order=args.order, ranks=args.ranks, gpu_batch_steps=args.batch, repeats=1,
                   precision='complex double', atol=1e-8, rtol=2e-6,
                   memory_note='Two-second whole-node samples, not per-process RSS or an exact peak',
                   cases=[], status='running')

    def save():
        path = root / 'summary.json'
        temporary = root / 'summary.json.tmp'
        temporary.write_text(json.dumps(summary, indent=2) + '\n')
        temporary.replace(path)

    save()
    try:
        # First validate the 32-rank shared-GPU launch on a small existing case.
        if args.order == 3:
            source = INPUTS / 'Diamond_NV_P1_10ppm_gcce/inputfiles/ccein_nstate0.json'
            pilot = root / 'mpi32-cce2-check'
            summary['current_case'] = '32-rank CCE2 MPI check'
            save()
            cpu = run_case(source, pilot / 'cpu', 2, 'cpu', args)
            gpu = run_case(source, pilot / 'cuda', 2, 'cuda', args)
            summary['mpi_check'] = dict(cpu=cpu, cuda=gpu, **compare(pilot / 'cpu', pilot / 'cuda', 1e-8, 2e-6))
            save()
        for case in sorted(INPUTS.glob('Diamond_*')):
            for state in (0, 1):
                source = case / 'inputfiles' / f'ccein_nstate{state}.json'
                name = f'{case.name}_nstate{state}'
                directory = root / name
                record = dict(case=name, backends={})
                summary['cases'].append(record)
                for backend in ('cpu', 'cuda', 'auto'):
                    summary.update(current_case=name, current_backend=backend)
                    save()
                    result = run_case(source, directory / backend, args.order, backend, args)
                    # Checking CPU against itself also rejects non-finite reference output.
                    error = compare(directory / 'cpu', directory / backend, 1e-8, 2e-6)
                    result.update(error)
                    if backend != 'cpu':
                        result['speedup'] = record['backends']['cpu']['wall_s'] / result['wall_s']
                    record['backends'][backend] = result
                    save()
                    print(f'order={args.order} {name} backend={backend} wall_s={result["wall_s"]:.3f} '
                          f'max_error={error["max_abs_error"]:.3g}', flush=True)
        if len(summary['cases']) != 8:
            raise ValueError('Expected exactly eight cases')
        summary['status'] = 'passed'
        save()
    except BaseException as error:
        summary.update(status='failed', error=f'{type(error).__name__}: {error}')
        save()
        raise


if __name__ == '__main__':
    main()
