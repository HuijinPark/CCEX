#!/usr/bin/env python3
"""Compare acceleration stages and MPI sizes with immutable inputs and raw checks."""
import argparse
import json
from pathlib import Path
import statistics
import sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[1]/'check_gpu'))   # benchmark.py, the helper the checks share
from benchmark import compare
from check_pcce import run_case
from check_acceleration import config_for
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'cluster_verification/check_pcce'))
from check_baths import resolve,sha,save


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cpu',type=Path,required=True);p.add_argument('--cuda',type=Path,required=True)
    p.add_argument('--before-cuda',type=Path,required=True);p.add_argument('--fixtures',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--rank-sweep',type=int,nargs='+',default=[1,2,4,8,16,32])
    args=p.parse_args();root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    fixture=args.fixtures.resolve();manifest=json.loads((fixture/'manifest.json').read_text())
    for name,want in manifest['fixture_sha256'].items():
        if sha(fixture/name)!=want:raise ValueError(f'Changed fixture {name}')
    cases={c['name']:c for c in manifest['cases']}
    selected=['diamond_13c_sK3_state0_pulse4','diamond_p1_15n_sK1_state0_pulse4',
        'diamond_p1_15n_sK2_state0_pulse1','diamond_p1_15n_sK2_state0_pulse4',
        'mos2_nuclear_sK2_state0_pulse1','mos2_nuclear_sK2_state0_pulse4',
        'synthetic_mixed_sK2_state0_pulse4','synthetic_mixed_sK3_state0_pulse4']
    binaries=dict(cpu=args.cpu.resolve(),cuda=args.cuda.resolve(),before_cuda=args.before_cuda.resolve())
    summary=dict(status='running',repeats=args.repeats,binary_sha256={k:sha(v) for k,v in binaries.items()},
        fixture_manifest_sha256=sha(fixture/'manifest.json'),stages=[],rank_sweep=[])
    def measure(name,config,ranks,variants,section):
        record=dict(case=name,ranks=ranks,config=config,trials=[])
        summary[section].append(record)
        for repeat in range(args.repeats):
            trial=dict(repeat=repeat,executions={},comparisons={});record['trials'].append(trial)
            order=variants[repeat%len(variants):]+variants[:repeat%len(variants)]
            directory=root/section/name/f'ranks{ranks}'/f'repeat{repeat}'
            for variant,binary,overrides in order:
                execution,log=run_case({**config,**overrides,'gpu_options':{**config.get('gpu_options',{}),**overrides.get('gpu_options',{})}},directory/variant,binary,ranks,600)
                if execution['returncode'] or not execution['profiles']:raise ValueError(f'{name}/{variant} failed')
                trial['executions'][variant]=execution
                save(root/'summary.json',summary)
            for variant,_,_ in variants:
                trial['comparisons'][variant]=compare(directory/'cpu',directory/variant,1e-8,2e-6)
                for key in ('cluster_listing_sha256','cluster_counts_by_dimension'):
                    if trial['executions'][variant][key]!=trial['executions']['cpu'][key]:raise ValueError(f'{name}/{variant}: {key} differs')
            print('PASS',section,name,'ranks',ranks,'repeat',repeat,
                {k:round(v['wall_s'],3) for k,v in trial['executions'].items()},flush=True)
        record['median_s']={v:statistics.median(t['executions'][v]['wall_s'] for t in record['trials']) for v,_,_ in variants}
        save(root/'summary.json',summary)
    for name in selected:
        c=resolve(cases[name]['config'],fixture);c['pcce_partition']='legacy'
        measure(name,c,2,[('cpu',binaries['cpu'],dict(backend='cpu')),
            ('before_cuda',binaries['before_cuda'],dict(backend='cuda')),
            ('cuda_sequential',binaries['cuda'],dict(backend='cuda',gpu_options=dict(evolution_batch_size=1))),
            ('cuda_batch',binaries['cuda'],dict(backend='cuda',gpu_options=dict(evolution_batch_size=8))),
            ('auto',binaries['cuda'],dict(backend='auto',gpu_options=dict(min_dim=256)))],'stages')
    sweep=[(name,dict(resolve(cases[name]['config'],fixture),pcce_partition='deterministic'))
        for name in ['diamond_13c_sK3_state0_pulse4','diamond_p1_15n_sK2_state0_pulse4']]
    c=config_for('Diamond_NV_natab_gcce',0);c['order']=4;sweep.append(('gcce4_13c',c))
    for name,c in sweep:
        reference=None
        for ranks in args.rank_sweep:
            measure(name,c,ranks,[('cpu',binaries['cpu'],dict(backend='cpu')),
                ('cuda_batch',binaries['cuda'],dict(backend='cuda',gpu_options=dict(evolution_batch_size=8))),
                ('auto',binaries['cuda'],dict(backend='auto',gpu_options=dict(min_dim=256)))],'rank_sweep')
            current=root/'rank_sweep'/name/f'ranks{ranks}'/'repeat0'/'cpu'
            if reference:compare(reference,current,1e-8,2e-6)
            else:reference=current
    summary['status']='passed';save(root/'summary.json',summary)


if __name__=='__main__':main()
