#!/usr/bin/env python3
"""The packed GPU batch (gpu_options.packed_max_dim) end to end without a card: its host side - the node shared memory,
the feeder thread, the windows - built with small_cluster_host_emulation.cpp in place of the card, on 3 ranks.

  python3 check_packed_emulation.py --output DIR      (CCEX_BUILD_DIR=<build> for a CMake build, CCEX_MPI=<mpirun>)

Builds DIR/main_emu.out from the objects of the current CPU build (test/unit_build_env.sh) with small_cluster_batch.cpp
compiled -DCCEX_SMALL_CLUSTER_EMULATE, then the regression cases of review 20261012-233058-452278:
  P3-01  a card that cannot be opened (CCEX_SMALL_EMULATE_OPEN_FAIL=1): backend auto runs every cluster on the CPU,
         byte-identical to packed_max_dim 0, whether the opening fails at once or after the ranks are done (400 ms) and
         with no clusters; an explicit packed_max_dim with backend cpu stops the run in all three;
  P3-04  the window keeps the 6 MiB budget at nstep 20000 (19 clusters, not 256), no window without clusters, and the
         batch off when one cluster's coherences alone pass the budget;
  and the batch on (GPU tolerance 1e-8 + 2e-6 |CPU| against the CPU).
"""
import argparse
import json
import os
from pathlib import Path
import re
import shlex
import subprocess

from benchmark import REPO, INPUTS, compare


def build(out):
    env_text = subprocess.run(['bash', str(REPO/'test/unit_build_env.sh')], check=True, capture_output=True, text=True).stdout
    unit = {k: v for k, v in re.findall(r"^(CCEX_UNIT_\w+)='(.*)'$", env_text, re.M)}
    env = dict(os.environ, **dict(p.split('=', 1) for p in shlex.split(unit['CCEX_UNIT_ENV'])))
    cxx = shlex.split(unit['CCEX_UNIT_CXX'])
    flags = shlex.split(unit['CCEX_UNIT_CXXFLAGS']) + shlex.split(unit['CCEX_UNIT_INCLUDES'])
    objects = [o for o in shlex.split(unit['CCEX_UNIT_OBJS']) if 'small_cluster_batch' not in Path(o).name]
    built = []
    log = (out/'build.log').open('w')        # the compiler's warnings
    for source, extra in [('src/calculation/gpu/small_cluster_batch.cpp', ['-DCCEX_SMALL_CLUSTER_EMULATE']),
                          ('test/gpu_verification/check_gpu/small_cluster_host_emulation.cpp', []),
                          ('src/main.cpp', [])]:
        target = out/(Path(source).stem+'.o')
        subprocess.run(cxx+flags+extra+['-c', str(REPO/source), '-o', str(target)], check=True, env=env,
                       stdout=log, stderr=subprocess.STDOUT)
        built.append(str(target))
    binary = out/'main_emu.out'
    subprocess.run(cxx+flags+built+objects+shlex.split(unit['CCEX_UNIT_LIBS'])+['-o', str(binary)], check=True, env=env,
                   stdout=log, stderr=subprocess.STDOUT)
    return binary


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--np', type=int, default=3)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    binary = build(out)
    case = INPUTS/'Diamond_NV_natab_gcce'
    original = json.loads((case/'inputfiles/ccein_nstate0.json').read_text())
    original.update(method='CCE', order=2, nstep=30, profile=True)
    original['gyrofile'] = str((case/original['gyrofile']).resolve())
    original['bathfile'] = [str((case/p).resolve()) for p in original['bathfile']]
    mpirun = os.environ.get('CCEX_MPI', 'mpirun')
    failures = []

    def run(name, overrides, env_extra=None, refused=None):
        directory = out/name
        directory.mkdir()
        config = dict(original, **overrides)
        config['outfile'] = str(directory/'coherence')
        (directory/'input.json').write_text(json.dumps(config))
        env = dict(os.environ, I_MPI_FABRICS='shm', **(env_extra or {}))
        with (directory/'run.log').open('w') as log:
            result = subprocess.run([mpirun, '-n', str(args.np), str(binary), '-f', str(directory/'input.json')],
                                    cwd=directory, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=600)
        text = (directory/'run.log').read_text()
        if refused:
            ok = result.returncode != 0 and refused in text
        else:
            ok = result.returncode == 0
        if not ok:
            failures.append(name)
            print(f'FAIL {name}: exit {result.returncode}\n{text[-1500:]}', flush=True)
        return ok, text, directory

    def same_bytes(a, b):
        files = sorted(p.name for p in a.iterdir() if p.name.startswith('coherence_'))
        return files and files == sorted(p.name for p in b.iterdir() if p.name.startswith('coherence_')) and \
            all((a/f).read_bytes() == (b/f).read_bytes() for f in files)

    def check(name, condition, what):
        print(('PASS ' if condition else 'FAIL ')+name+': '+what, flush=True)
        if not condition:
            failures.append(name)

    off = {'backend': 'auto', 'gpu_options': {'packed_max_dim': 0}}
    on = {'backend': 'auto', 'gpu_options': {'packed_max_dim': 4}}
    explicit_cpu = {'backend': 'cpu', 'gpu_options': {'packed_max_dim': 4}}
    no_card = {'CCEX_SMALL_EMULATE_OPEN_FAIL': '1'}
    _, _, reference = run('reference', off)

    # the batch on: GPU tolerance against the CPU
    ok, text, directory = run('packed-on', on)
    if ok:
        compare(reference, directory, 1e-8, 2e-6)
        check('packed-on', '[gpu-packed] on' in text, 'within 1e-8 + 2e-6 |CPU| of the CPU, batch on')

    # P3-01: auto without a card -> the CPU, byte for byte, at once, late and with no clusters; cpu + explicit -> refused
    for label, extra in [('fast', {}), ('slow', {'CCEX_SMALL_EMULATE_OPEN_MS': '400'})]:
        ok, text, directory = run('auto-no-card-'+label, on, dict(no_card, **extra))
        if ok:
            check('auto-no-card-'+label, same_bytes(reference, directory) and
                  '[gpu-packed] off: the GPU could not be opened' in text and 'clusters on the GPU 0,' in text,
                  'exit 0, every cluster on the CPU, curves byte-identical to packed_max_dim 0')
        ok, text, _ = run('cpu-no-card-'+label, explicit_cpu, dict(no_card, **extra), refused='Error: the packed GPU kernel')
        if ok:
            check('cpu-no-card-'+label, True, 'backend cpu with an explicit packed_max_dim and no card: refused')
    ok, text, _ = run('auto-no-card-zero', dict(on, rbath=1), no_card)
    if ok:
        check('auto-no-card-zero', '[gpu-packed] off: no clusters' in text and 'shared memory per rank' not in text,
              'no clusters: exit 0, no window')
    ok, text, _ = run('cpu-no-card-slow-zero-work', dict(explicit_cpu, order=1, rbath=8), dict(no_card, CCEX_SMALL_EMULATE_OPEN_MS='400'),
                      refused='Error: the packed GPU kernel')
    if ok:
        check('cpu-no-card-slow-zero-work', True, 'a few clusters done before the opening failed: still refused (after the join)')

    # P3-04: the window keeps the budget; none without clusters; off when one cluster passes it
    long = dict(on, order=1, rbath=12, nstep=20000, deltat=0.0001)
    _, _, long_reference = run('long-reference', dict(long, gpu_options={'packed_max_dim': 0}))
    ok, text, directory = run('long-window', long)
    if ok:
        compare(long_reference, directory, 1e-8, 2e-6)
        m = re.search(r'window (\d+), ([\d.]+) MiB shared memory per rank', text)
        check('long-window', m is not None and int(m.group(1)) == (6 << 20)//(16*20000) and float(m.group(2)) < 15,
              f'nstep 20000: window {m.group(1) if m else "?"}, {m.group(2) if m else "?"} MiB per rank, within the GPU tolerance')
    ok, text, _ = run('zero-clusters', dict(on, rbath=1, nstep=20000))
    if ok:
        check('zero-clusters', '[gpu-packed] off: no clusters' in text and 'shared memory per rank' not in text,
              'no clusters: no window, no shared memory')
    ok, text, directory = run('over-budget', dict(on, order=1, rbath=6, nstep=400000, deltat=0.00001))
    if ok:
        check('over-budget', 'pass the 6 MiB of a window' in text and 'shared memory per rank' not in text,
              'nstep 400000: the batch off, every cluster on the CPU')

    print(f'{"FAIL" if failures else "PASS"} check_packed_emulation: {len(failures)} failed {failures}')
    raise SystemExit(1 if failures else 0)


if __name__ == '__main__':
    main()
