#!/usr/bin/env python3
# HISTORICAL (2026-09-12..16): option names and binaries of that time; does not run with the current build.
# See README.md here; the current check is ../check_adaptive/verify.sh.
"""Check omitted group limits, old inputs, MPI, and bounded small-bath accuracy."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO/'test/gpu_verification/check_gpu'))
from check_pcce import run_case
from benchmark import compare, values


def write(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def check_plan(plan, cap, order):
    """Exhaustive dimension check, independent of the native graph walker."""
    groups = plan['groups']; dims = plan['site_dimensions']
    assert sorted(s for g in groups for s in g) == list(range(len(dims)))
    gd = [math.prod(dims[s] for s in g) for g in groups]
    q = plan['qubit_dimension']
    assert all(d*q <= cap for d in gd)
    edges = {(a, b) for a, b, _ in plan['edges']}
    checked = 0
    for size in range(1, min(order, len(groups))+1):
        for chosen in itertools.combinations(range(len(groups)), size):
            seen = {chosen[0]}
            while True:
                more = {s for s in chosen if any(tuple(sorted((s, t))) in edges for t in seen)}
                if more <= seen: break
                seen.update(more)
            if len(seen) == size:
                assert math.prod(gd[g] for g in chosen)*q <= cap, (chosen, gd, cap)
                checked += 1
    return dict(group_sites=[len(g) for g in groups], group_dimensions=gd,
                connected_combinations_checked=checked)


def curve_error(reference, candidate):
    files = sorted(reference.glob('*_wiDiv'))
    assert files and [f.name for f in files] == [f.name for f in sorted(candidate.glob('*_wiDiv'))]
    errors = []; scaled = []
    for f in files:
        before, after = values(f), values(candidate/f.name)
        assert len(before) == len(after)
        for a, b in zip(before, after):
            assert a[0] == b[0] and len(a) == len(b)
            for x, y in zip(a[1:], b[1:]):
                error = abs(x-y); errors.append(error); scaled.append(error/(1e-7+1e-5*abs(x)))
    assert errors
    return dict(max_abs_error=max(errors), max_scaled_error=max(scaled),
                points=len(errors), passes=max(scaled) <= 1)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, default=REPO/'bin/main.cpu.out')
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--backend', choices=('cpu', 'cuda'), default='cpu')
    p.add_argument('--skip-drivers', action='store_true')
    a = p.parse_args(); root = a.output.resolve(); root.mkdir(parents=True, exist_ok=False)
    base = json.loads((REPO/'test/cluster_verification/check_bundle/simple_ensemble.json').read_text())
    base.update(backend=a.backend, gpu_min_dim=1)
    report = dict(status='running', backend=a.backend, cases=[], accuracy=[], workflows=[],
                  binary_sha256=hashlib.sha256(a.binary.read_bytes()).hexdigest(),
                  baseline_sha256=hashlib.sha256(a.baseline.read_bytes()).hexdigest(),
                  limitations=['Small selected baths only; no full PL6/h-BN convergence claim.',
                               'Timings are single functional-test runs, not a performance benchmark.'])
    def save(): write(root/'summary.json', report)
    def run(name, cfg, old=False, reference=None, ranks=1, error=None):
        cfg = dict(cfg, backend='cpu' if old else a.backend)
        directory = root/name
        if cfg.get('cluster_mode') == 'adaptive':
            request = root/(name+'.request.json')
            write(request, dict(mode='run', response=str(directory/'response')))
            cfg['adaptive_work_file'] = str(request)
        execution, log = run_case(cfg, directory, (a.baseline if old else a.binary).resolve(), ranks, 300)
        row = dict(name=name, baseline=old, ranks=ranks, execution=execution)
        plan = None
        if error:
            assert execution['returncode'] != 0 and error in log, (name, log[-3000:])
            row['expected_error'] = error
        else:
            assert execution['returncode'] == 0, (name, log[-4000:])
            if cfg.get('cluster_mode') == 'adaptive':
                plan = json.loads((directory/'response.plan.json').read_text())
                row.update(check_plan(plan, cfg['adaptive_cluster_max_dim'], cfg['order']))
                if not old:
                    cap = cfg.get('adaptive_group_max_dim', cfg['adaptive_cluster_max_dim']//plan['qubit_dimension'])
                    mode = 'explicit' if 'adaptive_group_max_dim' in cfg else 'auto'
                    assert plan['dimension_limits']['group_max_dim'] == cap
                    assert plan['dimension_limits']['group_max_dim_source'] == mode
                    assert plan['merge_cutoff_rad_ms'] == cfg.get('adaptive_merge_cutoff', 0)
                    assert f'group_max_dim_source={mode}' in log
                    assert max(row['group_dimensions']) <= cap
                    assert all(int(d) <= cfg['adaptive_cluster_max_dim'] for d in execution['cluster_counts_by_dimension'])
            if reference: row['comparison'] = compare(reference, directory, 1e-8, 2e-6)
        report['cases'].append(row); save(); print('PASS', name, flush=True)
        return directory, plan

    def standard(cfg, n):
        return dict({k:v for k,v in cfg.items() if k != 'cluster_mode' and not k.startswith('adaptive_')},
                    order=n, rdip=1e6)

    try:
        # Omitted values match the pre-change program with explicit effective
        # limits. Check nonuniform groups and order-dependent boundary blocking.
        for method, order in itertools.product(('cce', 'gcce'), (2, 3, 4)):
            cfg = dict(base, method=method, order=order,
                       adaptive_cluster_max_dim=16 if method=='cce' else 49)
            explicit = dict(cfg, adaptive_group_max_dim=16, adaptive_merge_cutoff=0)
            name = f'limits-{method}-o{order}'
            ref, previous = run(name+'-before', explicit, old=True)
            current, plan = run(name, cfg, reference=ref)
            assert plan['groups'] == previous['groups'] and plan['edges'] == previous['edges']
            if order == 2:
                assert max(map(len, plan['groups'])) >= 3 and len(set(map(len, plan['groups']))) > 1
                _, mpi = run(name+'-mpi2', cfg, ranks=2, reference=current)
                assert mpi == plan
        # Explicit legacy inputs keep the same grouping and coherence.
        legacy = json.loads((REPO/'test/cluster_verification/check_bundle/ensemble.json').read_text())
        for method in ('cce', 'gcce'):
            cfg = dict(legacy, method=method, adaptive_cluster_max_dim=192)
            ref, before = run('legacy-'+method+'-before', cfg, old=True)
            _, after = run('legacy-'+method, cfg, reference=ref)
            assert before['groups'] == after['groups'] and before['edges'] == after['edges']
        # Full six-spin Hamiltonian is an independent oracle when the cap
        # permits the automatically formed single group. Include pulse/state modes.
        for method, state, pulses in itertools.product(('cce','gcce'), (0,1), (1,4)):
            cfg = dict(base, method=method, nstate=state, npulse=pulses,
                       adaptive_cluster_max_dim=64 if method=='cce' else 192)
            name = f'full-{method}-s{state}-p{pulses}'
            ref, _ = run(name+'-reference', standard(cfg, 6), old=True)
            current, plan = run(name, cfg)
            assert len(plan['groups']) == 1
            error = curve_error(ref, current); assert error['passes'], (name, error)
            report['accuracy'].append(dict(name=name, **error)); save()
        # PL6 six-spin region from the failure investigation, retaining its
        # field, pulse and full 129-point time grid. This is a selected bath.
        fixtures = REPO/'test/cluster_verification/check_bundle/pl6-convergence-20260913/fixtures'
        original = json.loads((fixtures/'pl6-original.json').read_text())
        cfg = dict(base, Qubit=original['Qubit'], bfield=57, rbath=40, rdip=10,
                   gyrofile=str(fixtures/'pl6-gyrofile'), nstep=129,
                   deltat=original['deltat'], seed=original['seed'])
        ids = [78,463,464,505,506,519]
        rows = (fixtures/'pl6-bath').read_text().splitlines()[1:]
        bath = root/'pl6-six.bath'; bath.write_text('7 0 0 0\n'+'\n'.join(rows[i] for i in ids)+'\n')
        cfg['bathfile'] = [str(bath)]
        report['pl6_sites'] = ids
        ref, _ = run('pl6-six-reference', standard(cfg, 6), old=True)
        for cap in (16,64):
            current, _ = run(f'pl6-six-cap{cap}', dict(cfg, adaptive_cluster_max_dim=cap))
            error = curve_error(ref,current)
            if cap == 64: assert error['passes'], error
            report['accuracy'].append(dict(name=f'pl6-six-cap{cap}', **error)); save()
        # Four real h-BN nuclei with mixed dimensions. These point-dipole tests
        # do not include the material's ab-initio hyperfine/quadrupole tensors.
        stored = REPO/'test/code_verification/bath_data/ensemble_data'
        xyz = list(map(float,(stored/'bath_h_BN_after_relax_Flat_1L_Defect').read_text().split()))
        rows = [line.split() for line in (stored/'bath_h_BN_after_relax_Flat_1L_1').read_text().splitlines()[1:]]
        nearest = sorted(enumerate(rows), key=lambda item:math.dist(xyz,list(map(float,item[1][:3]))))[:4]
        bath = root/'hbn-four.bath'; bath.write_text('5 0 0 0\n'+'\n'.join(' '.join(row) for _,row in nearest)+'\n')
        report['hbn_sites'] = dict(indices=[i for i,_ in nearest], species=[row[3] for _,row in nearest])
        for method in ('cce','gcce'):
            q = 3 if method=='gcce' else 1
            cfg = dict(base, method=method, bathfile=[str(bath)], gyrofile=str(stored/'Gyro_h_BN'),
                       Qubit=dict(nqubit=1, qubit=[dict(name='VB',spin=1,gyro=-17608.59705,xyz=xyz,alphams=1,betams=0)]),
                       deltat=.0001, npulse=4)
            ref,_ = run('hbn-'+method+'-reference', standard(cfg,4), old=True)
            for cap in (36,108):
                current,_ = run(f'hbn-{method}-cap{cap}', dict(cfg,adaptive_cluster_max_dim=cap*q))
                error = curve_error(ref,current)
                if cap == 108: assert error['passes'], error
                report['accuracy'].append(dict(name=f'hbn-{method}-cap{cap}', **error)); save()
        # Automatic *connections* still require their separate omission budget.
        for order in (2,3):
            cfg = dict(base,order=order,adaptive_connections='auto',adaptive_connection_budget=0)
            cfg.pop('rdip')
            ref,_ = run(f'auto-links-o{order}-before',dict(cfg,adaptive_group_max_dim=16,adaptive_merge_cutoff=0),old=True)
            run(f'auto-links-o{order}',cfg,reference=ref)
        for key,bad,message in [('adaptive_group_max_dim',0,'positive 32-bit'),
                                ('adaptive_group_max_dim',True,'finite JSON number'),
                                ('adaptive_group_max_dim',None,'finite JSON number'),
                                ('adaptive_merge_cutoff',-1,'nonnegative')]:
            run('reject-'+key+'-'+str(bad),dict(base,**{key:bad}),error=message)
        cfg = dict(base); cfg.pop('adaptive_cluster_max_dim')
        run('missing-cluster-cap',cfg,error='adaptive requires adaptive_cluster_max_dim')
        if not a.skip_drivers:
            # Exercise omitted keys in each Python workflow as well as native
            # propagation. Restrict validation explicitly; no expensive search.
            for policy in ('merge','error-gated','budgeted','auto-order'):
                inp = root/('driver-'+policy+'.json')
                write(inp,dict(base,adaptive_cluster_max_dim=64))
                dest = root/('driver-'+policy)
                command = [sys.executable,str(REPO/'test/cluster_verification/check_adaptive/oracle/adaptive_refine.py'),'--input',str(inp),
                           '--binary',str(a.binary.resolve()),'--output',str(dest),'--cwd',str(REPO),
                           '--policy',policy,'--atol','1e-7','--rtol','1e-5','--max-probes','64',
                           '--max-rounds','1','--validation-domains','2','--repair-rounds','0',
                           '--max-total-seconds','120','--max-order','3']
                with (root/('driver-'+policy+'.log')).open('w') as log:
                    result = subprocess.run(command,cwd=REPO,stdout=log,stderr=subprocess.STDOUT,timeout=180)
                summary = json.loads((dest/'summary.json').read_text())
                assert result.returncode==0, (policy,summary.get('status'),summary.get('error'))
                assert summary['status'] in ('full_bath_cluster','validated_full_bath_model')
                report['workflows'].append(dict(policy=policy,status=summary['status'])); save()
                print('PASS driver',policy,flush=True)
            # Both restricted real-bath inputs need reference domains larger
            # than their propagation cap. Never call an unperformed check a pass.
            for name in ('pl6-six-cap16','hbn-gcce-cap36'):
                cfg = json.loads((root/name/'input.json').read_text())
                cfg.pop('adaptive_work_file',None)
                inp = root/(name+'-audit.json'); write(inp,cfg)
                dest = root/(name+'-audit')
                command = [sys.executable,str(REPO/'test/cluster_verification/check_adaptive/oracle/adaptive_refine.py'),'--input',str(inp),
                           '--binary',str(a.binary.resolve()),'--output',str(dest),'--cwd',str(REPO),
                           '--policy','budgeted','--atol','1e-7','--rtol','1e-5','--max-probes','64',
                           '--max-rounds','0','--validation-domains','2','--repair-rounds','0','--timeout','60']
                with (root/(name+'-audit.log')).open('w') as log:
                    result = subprocess.run(command,cwd=REPO,stdout=log,stderr=subprocess.STDOUT,timeout=180)
                summary = json.loads((dest/'summary.json').read_text())
                assert result.returncode == 2 and summary['status'] == 'audit_budget_unverified'
                assert not summary['accuracy_certified']
                assert summary['cost']['probe_calls'] == 0
                assert summary['audit']['unique_domains'] <= 2
                assert summary['audit']['work_used'] <= summary['audit']['work_limit']
                report['workflows'].append(dict(name=name,policy='budgeted',status=summary['status'],
                                               accuracy_certified=False,audit=summary['audit']))
                save(); print('PASS bounded audit',name,flush=True)
        report['status'] = 'passed'; save()
    except BaseException as error:
        report.update(status='failed',error=f'{type(error).__name__}: {error}'); save(); raise


if __name__ == '__main__': main()
