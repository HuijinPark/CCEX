#!/usr/bin/env python3
"""How gpu_options.min_dim depends on the number of ranks that share one GPU (A3 of (internal notes)).

benchmark_min_dim.py measured one rank per GPU, and eight ranks on a problem too small to keep them busy. Here the
same P1 inputs (cut to the S bath sites nearest the NV) run at a fixed order K with R ranks on one GPU, S chosen as
the smallest number of sites with C(S, K) >= 4 R, so that at the deciding sizes every rank has at least four
clusters of the same size at the same time (CCEX computes one cluster size after the other on all ranks). cCCE:
K = 6 (dimensions 4 ... 4096); gCCE: K = 5 (3 ... 3072). Single sampled state, one pulse, 40 steps. backend=cpu and
backend=cuda (forced) from the same binary; per dimension the time per cluster summed over ranks, their ratio, and
the smallest dimension from which the GPU is faster at every larger one (only dimensions with at least 4 clusters
per rank count). A condition whose CPU and CUDA results differ beyond 1e-8 / 2e-6, or whose run fails, gets no crossover; the other conditions
still run, and the script exits 1 at the end.
usage: benchmark_gpu_sharing.py --binary BIN --mpirun MPIRUN --output DIR [--ranks 1 4 8 16 32] [--methods cce gcce]"""
import argparse
import json
from math import comb
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'check_gpu'))
from benchmark import compare
from benchmark_min_dim import prepare, run

ORDER = {'cce': 6, 'gcce': 5}


def sites_for(order, ranks):
    s = order
    while comb(s, order) < 4*ranks: s += 1
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--binary', type=Path, required=True); ap.add_argument('--mpirun', required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--ranks', type=int, nargs='+', default=[1, 4, 8, 16, 32])
    ap.add_argument('--methods', nargs='+', choices=sorted(ORDER), default=['cce', 'gcce'])
    ap.add_argument('--timeout', type=float, default=14400)
    args = ap.parse_args()
    binary = args.binary.resolve()
    if not binary.is_file(): sys.exit(f'no binary at {binary}')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    summary, invalid = [], []
    for method in args.methods:
        for ranks in args.ranks:
            order = ORDER[method]; sites = sites_for(order, ranks); name = f'{method}_r{ranks}'
            config = prepare(method, 1, sites, 1, 40, out/name/'input')
            config['order'] = order
            row = dict(name=name, method=method, ranks=ranks, sites=sites, order=order)
            try:
                for backend in ('cpu', 'cuda'):
                    row[backend+'_wall_s'], row[backend] = run(binary, config, out/name/backend, backend, ranks, args.mpirun, args.timeout)
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired, SystemExit) as e:
                # a failed run (e.g. the GPU memory of many ranks at once) is a result, not the end of the measurement
                row['failed'] = f'{backend}: {e}'; summary.append(row); invalid.append(name)
                print(f'{name}: sites {sites}, order {order}; {backend} run failed: {e} (see {out/name/backend}/run.log)', flush=True)
                (out/'summary.json').write_text(json.dumps(summary, indent=1, default=str)+'\n')
                continue
            try:
                row['difference'] = compare(out/name/'cpu', out/name/'cuda', 1e-8, 2e-6)
            except ValueError as e:
                row['difference'] = {'error': str(e)}; invalid.append(name)
            counts = {b: {d: v['clusters'] for d, v in row[b].items()} for b in ('cpu', 'cuda')}
            if counts['cpu'] != counts['cuda']:
                row['difference'] = dict(row.get('difference') or {}, work=f'cpu {counts["cpu"]} against cuda {counts["cuda"]}')
                if name not in invalid: invalid.append(name)
            ratio = {d: row['cpu'][d]['seconds']/row['cuda'][d]['seconds'] for d in row['cpu']
                     if d in row['cuda'] and row['cuda'][d]['seconds'] > 0}
            busy = [d for d in sorted(ratio) if row['cpu'][d]['clusters'] >= 4*ranks]
            row['speedup_by_dim'] = ratio; row['busy_dims'] = busy
            row['gpu_faster_from'] = None if name in invalid else \
                next((d for i, d in enumerate(busy) if all(ratio[e] > 1 for e in busy[i:])), None)
            summary.append(row)
            print(f'{name}: sites {sites}, order {order}; wall cpu {row["cpu_wall_s"]:.1f} s, cuda {row["cuda_wall_s"]:.1f} s; '
                  f'GPU faster from dim {row["gpu_faster_from"]}; difference {row["difference"]}', flush=True)
            for d in sorted(ratio):
                c, g = row['cpu'][d], row['cuda'][d]
                print(f'    dim {d:5d}  clusters {c["clusters"]:4d}{"" if d in busy else " (<4/rank)"}  per cluster cpu '
                      f'{c["seconds"]/c["clusters"]:.3e} s  cuda {g["seconds"]/g["clusters"]:.3e} s  cpu/cuda {ratio[d]:.2f}', flush=True)
            (out/'summary.json').write_text(json.dumps(summary, indent=1, default=str)+'\n')
    if invalid:
        sys.exit(f'CPU and CUDA results differ beyond the tolerance in: {", ".join(invalid)}')


if __name__ == '__main__':
    main()
