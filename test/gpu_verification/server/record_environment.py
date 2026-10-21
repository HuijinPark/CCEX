#!/usr/bin/env python3
"""Save the allocated server environment alongside the measured results."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--gpu', action='store_true')
    parser.add_argument('--binary', type=Path, nargs='+',
                        default=[Path('bin/main.out'), Path('bin/main_gpu.out')])
    args = parser.parse_args()
    record = dict(time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  hostname=platform.node(), python=platform.python_version(),
                  environment={k: os.environ[k] for k in (
                      'SLURM_JOB_ID', 'SLURM_JOB_NODELIST', 'SLURM_CPUS_PER_TASK',
                      'CUDA_VISIBLE_DEVICES', 'OMP_NUM_THREADS', 'LOADEDMODULES') if k in os.environ})
    commands = dict(cpu=['lscpu'], memory=['free', '-b'], compiler=['g++', '--version'],
                    cuda=['nvcc', '--version'], mpi=['mpirun', '--version'])
    if args.gpu:
        commands['gpu'] = ['nvidia-smi', '--query-gpu=index,name,memory.total,memory.used,driver_version,compute_cap',
                           '--format=csv']
    for key, command in commands.items():
        result = subprocess.run(command, capture_output=True, text=True, timeout=15, check=True)
        record[key] = (result.stdout + result.stderr).strip()
    record['binary_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in args.binary}
    with args.output.open('x') as stream:
        json.dump(record, stream, indent=2)
        stream.write('\n')


if __name__ == '__main__':
    main()
