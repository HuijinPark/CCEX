#!/usr/bin/env python3
"""The probe step across ranks and states (A1, 2026-10-02): probe i of state s runs on rank (i + s) % ranks, the ranks
do not wait for each other between states, and a node memory budget (probe_memory_gb) holds the large probes back.

    probe_rotation.py WORKDIR        # CCEX_BIN (default bin/main.out), CCEX_MPI (default mpirun)

The 14 13C nearest the NV of the natab fixture, three bath states, 24 probes of which three have 10 sites (dim 1024,
about 134 MB by the estimate, 8 dense matrices for a sampled cCCE state, so the budget counts them). Run on 1 rank, on 3, and on 3 with a 0.15 GB budget that
lets one large probe run at a time. Passes when the three give the same text for every (state, probe, step), every
probe once; on 3 ranks the large probes sit on the ranks the rotation says; and the small budget made a rank wait."""
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
BINARY = Path(os.environ.get('CCEX_BIN', REPO/'bin/main.out')).resolve()
MPIRUN = os.environ.get('CCEX_MPI') or 'mpirun'
FIXTURE = REPO/'test/code_verification/check_input_forms/inputs/Diamond_NV_natab_cce/inputfiles'
TEMPLATE = REPO/'test/cluster_verification/check_adaptive/reference/probe-step/input.json'
LARGE = [list(range(0, 10)), list(range(4, 14)), list(range(2, 12))]


def prepare(work):
    nv = [float(x) for x in (FIXTURE/'bath_Diamond_NV_natab_defect').read_text().split()[:3]]
    rows = [line.split() for line in (FIXTURE/'bath_Diamond_NV_natab_1').read_text().splitlines()[1:]]
    rows = sorted(rows, key=lambda r: math.dist([float(v) for v in r[:3]], nv))[:14]
    (work/'bath').write_text('15.0 0 0 0\n'+''.join('%.5f %.5f %.5f 13C\n' % tuple(float(r[i])-nv[i] for i in range(3))
                                                    for r in rows))
    probes = [[i] for i in range(14)]+[[i, i+1] for i in range(0, 13, 2)]+LARGE
    base = json.loads(TEMPLATE.read_text())
    base.update(gyrofile=str(FIXTURE/'gyro_13C'), bathfile=[str(work/'bath')], nstate=3, rbath=40, rdip=40,
                bundle=dict(base['bundle'], max_dim=2048), nstep=6, pulse={'count': 1}, profile=False, backend='cpu')
    return base, probes


def run(work, base, probes, name, ranks, budget):
    d = work/name
    d.mkdir()
    (d/'request.json').write_text(json.dumps(dict(mode='probe', response=str(d/'response'),
                                                  bundles=[[i] for i in range(14)], probes=probes)))
    cfg = dict(base, outfile=str(d/'coherence'), adaptive={'work_file': str(d/'request.json')})
    if budget:
        cfg['probe_memory_gb'] = budget
    (d/'input.json').write_text(json.dumps(cfg, indent=1))
    with open(d/'run.log', 'w') as log:
        rc = subprocess.run([MPIRUN, '-n', str(ranks), str(BINARY), '-f', str(d/'input.json')], cwd=REPO,
                            stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL).returncode
    if rc:
        sys.exit(f'{name}: exit {rc}, see {d}/run.log')
    values, owner = {}, {}
    for f in sorted(d.glob('response.rank-*.state-*.tsv')):
        m = re.search(r'rank-(\d+)\.state-(\d+)', f.name)
        assert m, f.name
        r, state = int(m.group(1)), int(m.group(2))
        for line in f.read_text().splitlines()[1:]:
            probe, step, real, imag = line.split('\t')
            key = (state, int(probe), int(step))
            if key in values:
                sys.exit(f'{name}: {key} written twice')
            values[key] = (real, imag)
            owner[(state, int(probe))] = r
    return values, owner, (d/'run.log').read_text()


def main(argv):
    if len(argv) != 1:
        print(__doc__); return 2
    work = Path(argv[0]).resolve()
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    base, probes = prepare(work)
    one, _, _ = run(work, base, probes, 'one-rank', 1, 0)
    three, owner, _ = run(work, base, probes, 'three-ranks', 3, 0)
    tight, _, log = run(work, base, probes, 'tight-budget', 3, 0.15)
    # probe -1 is the qubit alone, written by rank 0 (adaptiveProbe)
    expected = {(s, p, t) for s in (1, 2, 3) for p in range(-1, len(probes)) for t in range(base['nstep'])}
    failures = []
    if set(one) != expected:
        failures.append(f'one rank: {len(one)} values, {len(expected)} expected')
    if three != one:
        failures.append('three ranks differ from one')
    if tight != one:
        failures.append('three ranks with a 0.15 GB budget differ from one')
    first_large = len(probes)-len(LARGE)
    for s in (1, 2, 3):
        for p in range(first_large, len(probes)):
            if owner.get((s, p)) != (p+s) % 3:
                failures.append(f'state {s} probe {p} on rank {owner.get((s, p))}, not {(p+s) % 3}')
    if not re.search(r'\[probe-budget\] rank=\d+ waited', log):
        failures.append('the 0.15 GB budget made no rank wait')
    for line in failures:
        print('FAIL', line)
    if not failures:
        print(f'PASS probe rotation and memory budget: {len(one)} values alike on 1 rank, 3 ranks and 3 ranks with '
              f'a 0.15 GB budget; the large probes on rank (probe + state) % 3')
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
