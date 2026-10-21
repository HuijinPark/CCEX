#!/usr/bin/env python3
"""Prepare portable non-PL6 fixtures, then compare pCCE CPU/CUDA and old CPU.

This is a numerical regression suite, not a convergence study. Fixture settings
and all input changes are recorded. No simulation source is patched by this tool.
"""
import argparse
import collections
import copy
import hashlib
import json
import math
from pathlib import Path
import random
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'gpu_verification/check_gpu'))
from benchmark import REPO, compare, values, get, put
from check_pcce import run_case


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n')
    temporary.replace(path)


def bath_rows(path):
    return [line.split() for line in path.read_text().splitlines()[1:] if line.strip()]


def prepare(output):
    output.mkdir(parents=True, exist_ok=False)
    sources = {}

    def stage(source, target):
        source = source.resolve()
        shutil.copy2(source, output/target)
        sources[target] = dict(source=str(source.relative_to(REPO)), sha256=sha(source))
        return target

    families = []
    for name, folder, sizes, seed in (
            ('diamond_13c', 'Diamond_NV_natab_cce', [1, 2, 3], 73),
            ('diamond_p1_15n', 'Diamond_NV_P1_10ppm_cce', [1, 2], 97)):
        root = REPO/'test/code_verification/check_input_forms/inputs'/folder
        path = root/'inputfiles/ccein_nstate0.json'
        config = json.loads(path.read_text())
        for key in ('Qubit.file', 'configuration.gyrofile'):
            config = put(config, {key: stage(root/get(config, key), f'{name}-{key}')})
        config['configuration']['bathfile'] = [stage(root/config['configuration']['bathfile'][0], f'{name}-bath')]
        config['nstep'] = 17
        config['configuration']['seed'] = seed   # the configuration block (A64)
        families.append(dict(name=name, config=config, sizes=sizes,
            origin=str(path.relative_to(REPO)), input_sha256=sha(path),
            note='Original spatial bath, rbath, rdip, field, dt and defect tensors; method changed to pCCE, 17 time points, recorded seed and partition settings.'))

    root = REPO/'test/code_verification/bath_data/ensemble_data'
    config = dict(method='cCCE', cluster_mode='pcce', quantity='coherence', order=2, bfield=500,
        rbath=10, rdip=10, deltat=0.01, nstep=13)
    config['Qubit'] = dict(qubit=[dict(file=stage(root/'bath_MoS2_1L_Cimpurity_Defect', 'mos2-qubit'))], alpha_state=dict(ms=[1]), beta_state=dict(ms=[0]))
    config['configuration'] = dict(gyrofile=stage(root/'Gyro_MoS2', 'mos2-gyro'),
                                   bathfile=[stage(root/'bath_MoS2_1L_Cimpurity_1', 'mos2-bath')], seed=113)
    families.append(dict(name='mos2_nuclear', config=config, sizes=[1, 2],
        origin='test/code_verification/bath_data/NuclearSpinBath/cce.in',
        note='Stored 95Mo/97Mo/33S bath and legacy field/radii/dt; pCCE with 13 time points; point-dipole hyperfine, no quadrupole tensor in this fixture.'))

    # Reproducible artificial coordinates test different local Hilbert spaces
    # and a non-divisible site count; this is not a material model.
    gyro13 = (REPO/'test/code_verification/check_input_forms/inputs/Diamond_NV_natab_cce/inputfiles/gyro_13C').read_text()
    gyrobn = (root/'Gyro_h_BN').read_text().splitlines()
    (output/'mixed-gyro').write_text(gyro13.rstrip()+'\n'+'\n'.join(
        line for line in gyrobn if line.split()[0] in ('14N', '11B'))+'\n')
    rng = random.Random(20260911)
    points = []
    for species in ('13C', '14N', '11B', '13C', '14N', '11B', '13C'):
        while True:
            point = [round(rng.uniform(-10, 10), 5) for _ in range(3)]
            if math.dist(point, [0, 0, 0]) >= 4 and all(math.dist(point, old[:3]) >= 2 for old in points):
                break
        points.append([*point, species])
    (output/'mixed-bath').write_text('8 0 0 0\n'+''.join(
        ' '.join(map(str, row))+'\n' for row in points))
    (output/'mixed-qubit').write_text('0 0 0\n')
    config = dict(method='cCCE', cluster_mode='pcce', quantity='coherence', order=2,
        Qubit=dict(qubit=[dict(file='mixed-qubit')], alpha_state=dict(ms=[1]), beta_state=dict(ms=[0])), gyrofile='mixed-gyro', bathfile=['mixed-bath'],
        rbath=20, rdip=50, bfield=300,
        deltat=0.004, nstep=13, seed=211)
    families.append(dict(name='synthetic_mixed', config=config, sizes=[1, 2, 3],
        note='Artificial 7-site bath: spin 1/2, 1 and 3/2; seeded irregular geometry; no quadrupole. Gyros copied from stored 13C and h-BN tables.'))

    cases = []
    for family in families:
        config = family['config']
        origin = list(map(float, (output/config['Qubit']['qubit'][0]['file']).read_text().split()[:3]))
        inner, outer = config['rbath'] if isinstance(config['rbath'], list) else (0, config['rbath'])
        rows = [row for row in bath_rows(output/config['configuration']['bathfile'][0])
                if inner <= math.dist(list(map(float, row[:3])), origin) <= outer]
        spins = {line.split()[0]: float(line.split()[1]) for line in
                 (output/config['configuration']['gyrofile']).read_text().splitlines() if line.strip()}
        dimensions = {name: int(2*spin+1) for name, spin in spins.items()}
        for defect in config.get('Defect', {}).get('species', []):
            dimensions[defect['name']] *= math.prod(int(2*s+1) for s in defect.get('spins', []))
        family.update(loaded_sites=len(rows), species=dict(collections.Counter(row[3] for row in rows)),
                      site_dimensions=dimensions)
        variants = [(f'sK{size}', {'pcce.partition_size': size}, True) for size in family['sizes']]
        if family['name'] == 'synthetic_mixed':
            variants += [('one_center', {'pcce.partition_size': 7, 'pulse.count': 4}, False),
                         ('disconnected', {'pcce.partition_size': 2, 'rdip': 0.01, 'pulse.count': 4}, False),
                         ('empty_order4', {'pcce.partition_size': 2, 'order': 4, 'pulse.count': 4}, False),
                         ('random_init', {'pcce.partition_size': 2, 'pcce.kmeans_pp': False, 'pulse.count': 0, 'configuration.seed': 347}, True)]
        for label, overrides, baseline in variants:
            for state in (0, 1):
                for pulse in ([overrides['pulse.count']] if 'pulse.count' in overrides else [1, 4]):
                    c = put(config, {'method': 'cCCE', 'cluster_mode': 'pcce', 'order': 2, 'pcce.trials': 24,
                        'pcce.max_iterations': 40, 'pcce.kmeans_pp': True, 'pcce.print_trials': False,
                        'gpu_options.evolution_batch_steps': 5, 'gpu_options.min_dim': 16,
                        'profile': True, 'savemode': 'coherence', 'pulse.count': pulse, 'configuration.nstate': state})
                    c = put(c, overrides)
                    size = get(c, 'pcce.partition_size')
                    retained = rows[:len(rows)-len(rows)%size]
                    dims = sorted((dimensions[row[3]] for row in retained), reverse=True)
                    upper = math.prod(dims[:min(c['order']*size, len(dims))])
                    if upper > 1536:
                        raise ValueError(f'Fixture matrix dimension exceeds audit limit: {family["name"]} {label} {upper}')
                    # An ensemble physical reference checks cluster expansion,
                    # independently of old/new partition compatibility.
                    reference = None
                    if label == 'disconnected':
                        reference = dict(c, order=1)
                    elif state == 0 and (size == 1 or (family['name'] == 'synthetic_mixed' and c['order'] >= len(retained)//size)):
                        reference = put(c, {'cluster_mode': 'standard', 'pcce': None,
                                            'order': min(c['order'], len(retained)) if size == 1 else len(retained)})
                        if len(retained) != len(rows):
                            path = f'{family["name"]}-retained{len(retained)}'
                            (output/path).write_text(f'{len(retained)+1} 0 0 0\n'+''.join(' '.join(row)+'\n' for row in retained))
                            reference['configuration']['bathfile'] = [path]
                    cases.append(dict(name=f'{family["name"]}_{label}_state{state}_pulse{pulse}',
                        family=family['name'], config=c, compare_baseline=baseline, reference=reference,
                        loaded_sites=len(rows), retained_sites=len(retained), maximum_dimension_bound=upper))
    manifest = dict(families=families, cases=cases, source_files=sources,
        fixture_sha256={p.name: sha(p) for p in output.iterdir() if p.is_file()},
        scope='Non-PL6 pCCE numerical audit; no claim about physical convergence, T2 or full spatial ensemble averaging.')
    save(output/'manifest.json', manifest)
    print('Prepared', len(cases), 'cases:', [(f['name'], f['loaded_sites'], f['species']) for f in families])


def resolve(config, fixtures):
    config = copy.deepcopy(config)
    for key in ('Qubit.file', 'configuration.gyrofile'):
        config = put(config, {key: str(fixtures/get(config, key))})
    config['configuration']['bathfile'] = [str(fixtures/path) for path in config['configuration']['bathfile']]
    return config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--fixtures', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cpu', type=Path)
    parser.add_argument('--cuda', type=Path)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--ranks', type=int, default=2)
    parser.add_argument('--case', help='Optional case-name substring')
    parser.add_argument('--timeout', type=float, default=240)
    args = parser.parse_args()
    args.output = args.output.resolve()
    if args.prepare_only:
        prepare(args.output)
        return
    if not all((args.fixtures, args.cpu, args.cuda, args.baseline)):
        parser.error('--fixtures, --cpu, --cuda and --baseline are required')
    args.fixtures = args.fixtures.resolve()
    binaries = {name: getattr(args, name).resolve() for name in ('cpu', 'cuda', 'baseline')}
    manifest = json.loads((args.fixtures/'manifest.json').read_text())
    for name, expected in manifest['fixture_sha256'].items():
        if sha(args.fixtures/name) != expected:
            raise ValueError(f'Changed fixture: {name}')
    args.output.mkdir(parents=True, exist_ok=False)
    summary = dict(status='running', ranks=args.ranks, atol=1e-8, rtol=2e-6, repeats=1,
        families=manifest['families'], fixture_manifest_sha256=sha(args.fixtures/'manifest.json'),
        binary_sha256={name: sha(path) for name, path in binaries.items()}, cases=[])
    save(args.output/'summary.json', summary)
    for case in manifest['cases']:
        if args.case and args.case not in case['name']:
            continue
        record = dict(name=case['name'], family=case['family'], loaded_sites=case['loaded_sites'],
            retained_sites=case['retained_sites'], maximum_dimension_bound=case['maximum_dimension_bound'],
            executions={}, comparisons={}, status='running')
        summary['cases'].append(record)
        directory = args.output/case['name']
        try:
            config = resolve(case['config'], args.fixtures)
            for backend in (['baseline'] if case['compare_baseline'] else [])+['cpu', 'cuda']:
                result, log = run_case(dict(config, backend='cuda' if backend=='cuda' else 'cpu'),
                                      directory/backend, binaries[backend], args.ranks, args.timeout)
                record['executions'][backend] = result
                save(args.output/'summary.json', summary)
                if result['returncode'] or not result['profiles']:
                    raise RuntimeError(f'{backend} failed, returncode={result["returncode"]}')
                expected = 'cuda' if backend=='cuda' else 'cpu'
                if any(p['method']!='pcce' or p['backend']!=expected for p in result['profiles']):
                    raise ValueError('Unexpected method/backend profile')
                if max(int(p['dim']) for p in result['profiles']) > case['maximum_dimension_bound']:
                    raise ValueError('Unexpected matrix dimension')
                compare(directory/backend, directory/backend, 1e-8, 2e-6)
            for before, after in [('baseline', 'cpu'), ('cpu', 'cuda')]:
                if before not in record['executions']:
                    continue
                a, b = record['executions'][before], record['executions'][after]
                if a['cluster_listing_sha256'] != b['cluster_listing_sha256'] or a['cluster_counts_by_dimension'] != b['cluster_counts_by_dimension']:
                    raise ValueError(f'{before}/{after} cluster partition changed')
                record['comparisons'][f'{before}_{after}'] = compare(directory/before, directory/after, 1e-8, 2e-6)
            if case['reference']:
                result, _ = run_case(dict(resolve(case['reference'], args.fixtures), backend='cpu'),
                                     directory/'physical_reference', binaries['cpu'], args.ranks, args.timeout)
                record['executions']['physical_reference'] = result
                if result['returncode']:
                    raise RuntimeError('Physical reference failed')
                record['comparisons']['physical_reference'] = compare(directory/'physical_reference', directory/'cpu', 1e-8, 2e-6)
            files = list((directory/'cpu').glob('*_wiDiv'))
            record['signal_change_from_initial'] = max(abs(row[1]-rows[0][1]) for path in files
                for rows in [values(path)] for row in rows)
            record['status'] = 'passed'
            print('PASS',case['name'], 'dim<=',case['maximum_dimension_bound'],
                  'cpu_s=',round(record['executions']['cpu']['wall_s'],3),
                  'cuda_s=',round(record['executions']['cuda']['wall_s'],3), flush=True)
        except Exception as error:
            record.update(status='failed', error=f'{type(error).__name__}: {error}')
            print('FAIL',case['name'],record['error'],flush=True)
        save(args.output/'summary.json', summary)
    summary['status'] = 'passed' if summary['cases'] and all(c['status']=='passed' for c in summary['cases']) else 'failed'
    save(args.output/'summary.json', summary)
    if summary['status'] != 'passed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
