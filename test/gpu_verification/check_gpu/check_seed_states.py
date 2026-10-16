#!/usr/bin/env python3
"""Exercise seed sensitivity, state replay and CPU/CUDA/MPI consistency in main."""
import argparse
import ctypes
import hashlib
import json
from pathlib import Path
import re

from benchmark import REPO, compare, put
from check_pcce import run_case


def sampled(log):
    def values(label, text):
        return [(int(i), tuple(float(x) for x in value.strip().strip('[]').split()))
                for i, value in re.findall(
                    r'^\s*'+label+r'\[(\d+)\](?:\.state)?\s*:\s*([^\n]+)', text, re.M)]
    def state_blocks(kind, label):
        sections = re.findall(r'>> Randomize '+kind+r'-bath states\.\.\.(.*?)(?=\n\s*>>|\Z)', log, re.S)
        return [item for section in sections for item in values(label, section)]
    # P1 verbose logs also print each propagated cluster's renumbered spins.
    # Only inspect the original main/subbath sampling blocks, not those subsets.
    return dict(bath=state_blocks('main', 'bath'), subbath=state_blocks('sub', 'subbath'),
                axes=values('paxes', log))


def random_reference(seed, size, states, defect_axes=0):
    # The executable uses this platform's C rand(), not Python's RNG. Generate
    # an independent expected stream, including defect-axis draws before states.
    libc = ctypes.CDLL(None)
    libc.srand.argtypes = [ctypes.c_uint]
    libc.rand.restype = ctypes.c_int
    libc.srand(seed)
    axes = [(i, (float(libc.rand() % defect_axes + 1),)) for i in range(size)] if defect_axes else []
    bath, subbath = [], []
    for _ in range(states):
        bath.extend((i, (0.5 - libc.rand() % 2,)) for i in range(size))
        if defect_axes:
            subbath.extend((i, (0.5 - libc.rand() % 2,)) for i in range(size))
    return dict(bath=bath, subbath=subbath, axes=axes)


def check_stream(actual, expected):
    for key in ('bath', 'subbath'):
        if not expected[key]:
            assert not actual[key], (key, actual[key])
            continue
        count = len(expected[key])
        assert actual[key] and len(actual[key]) % count == 0, (key, len(actual[key]), count)
        # An rdip-omitted adaptive audit may print the states once during its
        # checks and again when those same states are restored for the final run.
        for i in range(0, len(actual[key]), count):
            assert actual[key][i:i+count] == expected[key], (key, actual[key], expected[key])
    assert actual['axes'] == expected['axes'], (actual['axes'], expected['axes'])


def main(args):
    args.output.mkdir(parents=True, exist_ok=False)
    base = json.loads((REPO/'test/gpu_verification/check_gpu/examples/pcce_ensemble.json').read_text())
    base = put(base, {'backend': 'cpu', 'configuration.nstate': 3, 'nstep': 5, 'pcce.trials': 5, 'pcce.print_trials': False})
    # the same input without the pcce partition: cluster_mode standard, method cCCE or gCCE (A47)
    standard = lambda method: put(base, {'method': method, 'cluster_mode': 'standard', 'pcce': None})
    report = dict(status='running', binary_sha256={}, cases=[], negative_control=None)
    for name in ('cpu', 'cuda', 'baseline'):
        binary = getattr(args, name)
        if binary: report['binary_sha256'][name] = hashlib.sha256(binary.read_bytes()).hexdigest()

    def save():
        (args.output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')

    def run(label, cfg, binary=None, ranks=1):
        directory = args.output/label
        execution, log = run_case(cfg, directory, binary or args.cpu, ranks)
        assert execution['returncode'] == 0, (label, log[-3000:])
        compare(directory, directory, 1e-8, 2e-6)  # Nonempty, finite coherence outputs.
        states = sampled(log)
        assert states['bath'], label
        return directory, states, log

    try:
        if args.baseline:
            old = [run(f'before/seed-{s}', put(standard('gCCE'), {'configuration.seed': s}), args.baseline)[1]
                   for s in (19, 7, 12345)]
            assert old[0] == old[1] == old[2], 'Baseline no longer reproduces the reported seed bug'
            check_stream(old[0], random_reference(0, 6, 3))
            report['negative_control'] = 'Before binary ignores all three seeds and follows srand(0)'
        p1 = json.loads((REPO/'test/gpu_verification/check_gpu/examples/p1_gcce4_single.json').read_text())
        p1.update(backend='cpu', order=1, nstep=4, profile=True)
        p1['configuration']['nstate'] = 3
        cases = [('cce', standard('cCCE'), 6, False),
                 ('gcce', standard('gCCE'), 6, False),
                 ('pcce', put(base, {'pcce.partition': 'deterministic'}), 6, False),
                 ('pcce-legacy', put(base, {'pcce.partition': 'legacy'}), 6, False),
                 ('adaptive', put(standard('cCCE'), {'cluster_mode': 'bundle', 'bundle.max_dim': 16}), 6, False),
                 ('adaptive-checked', put(standard('gCCE'), {'cluster_mode': 'bundle', 'bundle.max_dim': 48,
                                                             'deltat': 0.005, 'rdip': None}), 6, False),
                 ('p1-gcce', p1, 4, True),
                 ('p1-cce', put(p1, {'method': 'cCCE'}), 4, True)]
        for name, cfg, size, defect in cases:
            observed = []
            row = dict(name=name, checks=[])
            for seed in (19, 7, 12345):
                current = put(cfg, {'configuration.seed': seed})
                directory, states, _ = run(f'{name}/seed-{seed}', current)
                if name != 'pcce-legacy':
                    axes = cfg['Defect']['species'][0]['n_axis_states'] if defect else 0
                    check_stream(states, random_reference(seed, size, 3, axes))
                observed.append(states)
                row['checks'].append(dict(seed=seed, sampled_states=states))
                if seed == 19: reference, first = directory, states
            assert len({json.dumps(s, sort_keys=True) for s in observed}) == 3, name
            # Inspect random defect axes and subspin states separately as well.
            if defect:
                for key in ('axes', 'subbath'):
                    assert len({json.dumps(s[key]) for s in observed}) > 1, (name, key)
            # Legacy consumes rank-dependent rand() calls while partitioning;
            # only deterministic pCCE promises MPI-independent sampled states.
            variants = [('repeat', args.cpu, 1, 'cpu')]
            if name != 'pcce-legacy': variants.append(('mpi4', args.cpu, 4, 'cpu'))
            if args.cuda:
                variants.extend([('cuda', args.cuda, 1, 'cuda'), ('auto', args.cuda, 1, 'auto')])
            for label, binary, ranks, backend in variants:
                directory, states, _ = run(f'{name}/{label}', put(cfg, {'configuration.seed': 19, 'backend': backend}), binary, ranks)
                assert states == first, (name, label, 'state replay')
                row['checks'].append(dict(label=label, ranks=ranks,
                                         comparison=compare(reference, directory, 1e-8, 2e-6)))
            if args.baseline:
                # Historical-output compatibility is checked using the actual
                # old C RNG seed, separately from the changed explicit-seed semantics.
                before, _, _ = run(f'{name}/zero-before', put(cfg, {'configuration.seed': 0}), args.baseline)
                after, _, _ = run(f'{name}/zero-after', put(cfg, {'configuration.seed': 0}))
                row['zero_seed_compatibility'] = compare(before, after, 1e-8, 2e-6)
            report['cases'].append(row); save()
            print('PASS', name, 'seed sensitivity, expected states, replay and supported backends/ranks', flush=True)
        # Missing-seed replay must also replay actual sampled states.
        cfg = put(standard('gCCE'), {'configuration.seed': None})
        reference, states, log = run('generated/initial', cfg, ranks=2)
        seeds = re.findall(r'\[seed\] generated on rank 0 and shared with all ranks: (\d+),', log)
        assert len(seeds) == 1, seeds
        seed = int(seeds[0]); check_stream(states, random_reference(seed, 6, 3))
        report['generated_seed'] = seed
        for ranks in (1, 4):
            directory, replay, _ = run(f'generated/replay-{ranks}', put(cfg, {'configuration.seed': seed}), ranks=ranks)
            assert replay == states
            compare(reference, directory, 1e-8, 2e-6)
        report['status'] = 'passed'; save()
        print('PASS generated seed replays sampled states across MPI ranks', flush=True)
    except BaseException as error:
        report.update(status='failed', error=f'{type(error).__name__}: {error}'); save(); raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpu', type=Path, required=True)
    parser.add_argument('--cuda', type=Path)
    parser.add_argument('--baseline', type=Path, help='Optional pre-fix binary for the negative control')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    for key, value in vars(args).items():
        if value is not None: setattr(args, key, value.resolve())
    main(args)
