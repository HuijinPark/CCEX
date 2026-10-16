#!/usr/bin/env python3
"""Integration regression for gCCE4 accuracy, cluster batching and lazy auto."""
import argparse
import json
from pathlib import Path
import re
from benchmark import INPUTS,compare,get,put
from check_pcce import run_case


def config_for(family,state):
    root=INPUTS/family
    c=json.loads((root/'inputfiles'/f'ccein_nstate{state}.json').read_text())
    for key in ('configuration.gyrofile','Qubit.file'):
        if get(c,key):c=put(c,{key:str(root/get(c,key))})
    c['configuration']['bathfile']=[str(root/p) for p in c['configuration']['bathfile']]
    c=put(c,{'profile':True,'savemode':'coherence','propagator':'eigen','gpu_options':dict(get(c,'gpu_options',{}),evolution_batch_steps=4)})
    return c


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cpu',type=Path,required=True);p.add_argument('--cuda',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--ranks',type=int,default=2)
    p.add_argument('--rank-reference',type=int,default=1)
    args=p.parse_args();root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    records=[]
    cases=[]
    for state in (0,1):
        c=config_for('Diamond_NV_natab_gcce',state);c.update(order=4)
        cases.append((f'gcce4_13c_state{state}',c))
    c=put(config_for('Diamond_NV_natab_gcce',1),{'order':2,'evolution':'matrix','pulse.count':4,'nstep':9})
    cases.append(('gcce_single_matrix_pulse4',c))
    for method in ('cce','gcce'):
        for state in (0,1):
            c=put(config_for('Diamond_NV_P1_10ppm_'+method,state),{'order':2,'nstep':9,'pulse.count':4})
            cases.append((f'p1_{method}_state{state}',c))
    for name,c in cases:
        record=dict(case=name,executions={},comparisons={})
        for variant,overrides in [('cpu',dict(backend='cpu')),('cuda_sequential',dict(backend='cuda',gpu_options=dict(get(c,'gpu_options'),evolution_batch_size=1))),
                ('cuda_batch',dict(backend='cuda',gpu_options=dict(get(c,'gpu_options'),evolution_batch_size=8))),('auto',dict(backend='auto'))]:
            binary=args.cpu if variant=='cpu' else args.cuda
            execution,log=run_case(dict(c,**overrides),root/name/variant,binary.resolve(),args.ranks,300)
            if execution['returncode'] or not execution['profiles']:raise ValueError(f'{name}/{variant} failed')
            error=compare(root/name/'cpu',root/name/variant,1e-8,2e-6)
            if variant!='cpu':
                for key in ('cluster_listing_sha256','cluster_counts_by_dimension'):
                    if execution[key]!=record['executions']['cpu'][key]:raise ValueError(f'{name}/{variant}: {key} differs')
            record['executions'][variant]=execution;record['comparisons'][variant]=error
            if variant=='cuda_batch' and '[gpu-batch]' not in log:raise ValueError(f'Batch path not exercised: {name}')
            if variant=='auto' and max(map(int,execution['cluster_counts_by_dimension']))<get(c,'gpu_options').get('min_dim',768):
                if re.search(r'\[backend\].*CUDA device=',log) or '[gpu-memory]' in log:
                    raise ValueError('Small auto run unnecessarily initialized CUDA')
            print('PASS',name,variant,error,flush=True)
        if name.startswith('gcce4_') and args.ranks!=args.rank_reference:
            directory=root/name/f'cpu_ranks{args.rank_reference}'
            execution,_=run_case(dict(c,backend='cpu'),directory,args.cpu.resolve(),args.rank_reference,300)
            if execution['returncode']:raise ValueError('gCCE rank reference failed')
            record['rank_reference']=dict(ranks=args.rank_reference,execution=execution,
                comparison=compare(root/name/'cpu',directory,1e-8,2e-6))
            print('PASS',name,'MPI rank consistency',record['rank_reference']['comparison'],flush=True)
        records.append(record)
        (root/'summary.json').write_text(json.dumps(dict(status='running',ranks=args.ranks,cases=records),indent=2)+'\n')
    (root/'summary.json').write_text(json.dumps(dict(status='passed',ranks=args.ranks,cases=records),indent=2)+'\n')


if __name__=='__main__':main()
