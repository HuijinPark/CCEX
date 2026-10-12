#!/usr/bin/env python3
"""Check pCCE global trial totals, rank-independent partitions and bath states."""
import argparse
import json
from pathlib import Path
import re
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'gpu_verification/check_gpu'))
from benchmark import REPO,compare
from check_pcce import run_case
from partition_log import check_search


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpu',type=Path,required=True)
    parser.add_argument('--cuda',type=Path)
    parser.add_argument('--baseline',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--ranks',type=int,nargs='+',default=[1,2,3,8])
    args=parser.parse_args()
    root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    config=json.loads((REPO/'test/gpu_verification/check_gpu/examples/pcce_ensemble.json').read_text())
    config['gyrofile']=str(REPO/config['gyrofile'])
    config['bathfile']=[str(REPO/p) for p in config['bathfile']]
    config.update(iter_detail=True,backend='cpu')
    records=[]
    for size,state,trials,seed,plus in [(1,0,5,19,True),(2,2,5,19,True),
            (2,2,17,73,False),(3,2,1,97,True),(2,0,1,0,False)]:
        case=f'sK{size}_state{state}_trials{trials}_seed{seed}_plus{int(plus)}'
        previous=None;reference=None;search_reference=None
        for ranks in args.ranks:
            directory=root/case/f'ranks{ranks}'
            run,log=run_case(dict(config,sK=size,nstate=state,max_trial=trials,seed=seed,kmeans_pp=plus),directory,args.cpu.resolve(),ranks,120)
            if run['returncode']:raise RuntimeError(f'Execution failed: {case} ranks={ranks}')
            audit=check_search(log,trials,ranks)
            if search_reference is None:search_reference=audit
            elif any(audit[key]!=search_reference[key] for key in ('best_trial','search')):
                raise ValueError(f'Rank-dependent search work or winning trial: {case}')
            if previous:
                for key in ('cluster_listing_sha256','cluster_counts_by_dimension'):
                    if run[key]!=previous[key]:raise ValueError(f'Rank-dependent {key}: {case}')
                error=compare(reference,directory,1e-8,2e-6)
            else: reference=directory;error=compare(directory,directory,1e-8,2e-6)
            previous=run
            records.append(dict(case=case,ranks=ranks,execution=run,comparison=error,trial_audit=audit))
            print('PASS',case,'ranks',ranks,flush=True)
    if args.baseline:
        for ranks in [1,2]:
            for mode,binary in [('before',args.baseline),('legacy',args.cpu)]:
                # Old main reset the C RNG to zero regardless of the input seed.
                # Seed 0 compares the same historical random stream explicitly.
                result,_=run_case(dict(config,pcce_partition='legacy',nstate=2,seed=0),root/f'legacy{ranks}'/mode,binary.resolve(),ranks,120)
                if result['returncode']:raise RuntimeError('Legacy execution failed')
            error=compare(root/f'legacy{ranks}'/'before',root/f'legacy{ranks}'/'legacy',1e-8,2e-6)
            records.append(dict(case='legacy_compatibility',ranks=ranks,comparison=error))
            print('PASS legacy compatibility',ranks,flush=True)
    if args.cuda:
        # Only one rank receives the single six-spin cluster. Lazy initialization
        # must not wait for ranks that never encounter a GPU-eligible cluster.
        c=dict(config,sK=3,max_trial=1,nstate=2,gpu_options=dict(config.get('gpu_options',{}),min_dim=64))
        ranks=max(args.ranks)
        cpu,_=run_case(dict(c,backend='cpu'),root/'uneven_auto/cpu',args.cpu.resolve(),ranks,120)
        gpu,log=run_case(dict(c,backend='auto'),root/'uneven_auto/auto',args.cuda.resolve(),ranks,120)
        if cpu['returncode'] or gpu['returncode']:raise RuntimeError('Uneven auto execution failed')
        devices=re.findall(r'\[backend\] rank=(\d+) CUDA device=',log)
        if len(devices)!=1:raise ValueError(f'Expected exactly one lazy GPU initialization: {devices}')
        error=compare(root/'uneven_auto/cpu',root/'uneven_auto/auto',1e-8,2e-6)
        records.append(dict(case='uneven_auto',ranks=ranks,gpu_ranks=devices,comparison=error))
        print('PASS uneven lazy initialization',ranks,flush=True)
    (root/'summary.json').write_text(json.dumps(dict(status='passed',records=records),indent=2)+'\n')


if __name__=='__main__':main()
