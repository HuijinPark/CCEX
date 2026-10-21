#!/usr/bin/env python3
"""Compare CPU/CUDA/auto MPI scaling for one immutable pCCE bath fixture."""
import argparse
import json
from pathlib import Path
import statistics
import sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[1]/'check_gpu'))   # benchmark.py, the helper the checks share
from benchmark import REPO,compare,put
from check_pcce import run_case
sys.path.insert(0,str(REPO/'test/cluster_verification/check_pcce'))
from check_baths import sha,save,resolve


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fixtures',type=Path,required=True)
    p.add_argument('--case',default='mos2_nuclear_sK2_state0_pulse1')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ranks',type=int,nargs='+',default=[1,2,4,8,16,32])
    p.add_argument('--repeats',type=int,default=3)
    args=p.parse_args();root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    fixture=args.fixtures.resolve();manifest=json.loads((fixture/'manifest.json').read_text())
    for name,want in manifest['fixture_sha256'].items():
        if sha(fixture/name)!=want:raise ValueError(f'Fixture changed: {name}')
    case=next(c for c in manifest['cases'] if c['name']==args.case)
    config=resolve(case['config'],fixture)
    config=put(config,{'pcce.partition':'deterministic','gpu_options.min_dim':256})
    binaries=dict(cpu=REPO/'bin/main.out',cuda=REPO/'bin/main_gpu.out')
    summary=dict(status='running',case=args.case,config=config,repeats=args.repeats,rank_sweep=[],
        fixture_manifest_sha256=sha(fixture/'manifest.json'),binary_sha256={k:sha(v) for k,v in binaries.items()})
    reference=None;reference_execution=None
    try:
        for ranks in args.ranks:
            record=dict(ranks=ranks,trials=[]);summary['rank_sweep'].append(record)
            for repeat in range(args.repeats):
                variants=['cpu','cuda','auto'];variants=variants[repeat%len(variants):]+variants[:repeat%len(variants)]
                directory=root/f'ranks{ranks}'/f'repeat{repeat}'
                trial=dict(repeat=repeat,executions={},comparisons={});record['trials'].append(trial)
                for variant in variants:
                    summary.update(current_ranks=ranks,current_repeat=repeat,current_backend=variant);save(root/'summary.json',summary)
                    result,_=run_case(dict(config,backend=variant),directory/variant,
                        binaries['cpu' if variant=='cpu' else 'cuda'],ranks,600)
                    if result['returncode'] or not result['profiles']:raise RuntimeError(f'{ranks}/{variant} failed')
                    trial['executions'][variant]=result;save(root/'summary.json',summary)
                for variant in variants:
                    trial['comparisons'][variant]=compare(directory/'cpu',directory/variant,1e-8,2e-6)
                    for key in ('cluster_listing_sha256','cluster_counts_by_dimension'):
                        if trial['executions'][variant][key]!=trial['executions']['cpu'][key]:raise ValueError(f'{key} changed')
                print('PASS',args.case,'ranks',ranks,'repeat',repeat,{k:round(v['wall_s'],3) for k,v in trial['executions'].items()},flush=True)
            current=root/f'ranks{ranks}'/'repeat0/cpu'
            if reference:
                record['rank_comparison']=compare(reference,current,1e-8,2e-6)
                for key in ('cluster_listing_sha256','cluster_counts_by_dimension'):
                    if record['trials'][0]['executions']['cpu'][key]!=reference_execution[key]:raise ValueError(f'Rank-dependent {key}')
            else:reference=current;reference_execution=record['trials'][0]['executions']['cpu']
            record['median_s']={v:statistics.median(t['executions'][v]['wall_s'] for t in record['trials']) for v in ['cpu','cuda','auto']}
            save(root/'summary.json',summary)
        summary['status']='passed';save(root/'summary.json',summary)
    except BaseException as error:
        summary.update(status='failed',error=f'{type(error).__name__}: {error}');save(root/'summary.json',summary);raise


if __name__=='__main__':main()
