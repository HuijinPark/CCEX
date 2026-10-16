#!/usr/bin/env python3
"""Repeated PL6 comparisons using unchanged inputs from the previous server run."""
import argparse
import json
from pathlib import Path
import statistics
from benchmark import REPO,compare,get,put
from check_pcce import run_case
import sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[1]/'check_gpu'))   # benchmark.py, the helper the checks share
sys.path.insert(0,str(REPO/'test/cluster_verification/check_pcce'))
from check_baths import sha,save


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-run',type=Path,required=True)
    parser.add_argument('--before-cuda',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--ranks',type=int,default=32)
    parser.add_argument('--repeats',type=int,default=3)
    parser.add_argument('--cases',nargs='+',default=['ccein_fine_N1_sK1_rd6_longwindow','ccein_fine_N1_sK3_rd10'])
    args=parser.parse_args();root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    binaries=dict(cpu=REPO/'bin/main.out',cuda=REPO/'bin/main_gpu.out',before_cuda=args.before_cuda.resolve())
    summary=dict(status='running',ranks=args.ranks,repeats=args.repeats,
        binary_sha256={k:sha(v) for k,v in binaries.items()},cases=[],
        note='Physical settings, seed and time grid from previous PL6 config-1 inputs. Explicit legacy partition preserves the old 32-rank partition. Compare raw outputs; timings include fresh MPI/CUDA startup.')
    variants=[('cpu','cpu',dict(backend='cpu')),('before_cuda','before_cuda',dict(backend='cuda')),
        ('before_auto','before_cuda',dict(backend='auto',gpu_options=dict(min_dim=256))),
        ('cuda_sequential','cuda',dict(backend='cuda',gpu_options=dict(evolution_batch_size=1))),
        ('cuda_batch','cuda',dict(backend='cuda',gpu_options=dict(evolution_batch_size=8))),
        ('auto','cuda',dict(backend='auto',gpu_options=dict(min_dim=256)))]
    try:
        for name in args.cases:
            source=args.baseline_run.resolve()/name/'cpu/input.json'
            config=put(json.loads(source.read_text()),{'pcce.partition':'legacy'})
            # the configuration block (A64): a baseline input written before it is put through
            # example/tools/InputKeys/convert_input_keys.py first
            assert get(config,'cluster_mode')=='pcce' and get(config,'configuration.nstate',0)==0 and config['order']==2
            data=[get(config,'configuration.gyrofile'),*get(config,'configuration.bathfile')]
            if get(config,'Qubit.file'):data.append(get(config,'Qubit.file'))
            record=dict(case=name,source_input=str(source),source_input_sha256=sha(source),
                input_data_sha256={p:sha(Path(p)) for p in data},config=config,trials=[])
            summary['cases'].append(record)
            for repeat in range(args.repeats):
                trial=dict(repeat=repeat,executions={},comparisons={});record['trials'].append(trial)
                directory=root/name/f'repeat{repeat}'
                order=variants[repeat%len(variants):]+variants[:repeat%len(variants)]
                for variant,binary,overrides in order:
                    summary.update(current_case=name,current_repeat=repeat,current_variant=variant);save(root/'summary.json',summary)
                    result,log=run_case({**config,**overrides,'gpu_options':{**get(config,'gpu_options',{}),**overrides.get('gpu_options',{})}},directory/variant,binaries[binary],args.ranks,600)
                    if result['returncode'] or not result['profiles']:raise RuntimeError(f'{name}/{variant} failed')
                    if variant=='cuda_batch' and '[gpu-batch]' not in log:raise RuntimeError('Cluster batching was not exercised')
                    if variant=='auto' and ('[gpu-memory]' in log or '[backend] rank=' in log):raise RuntimeError('Small PL6 auto unnecessarily initialized CUDA')
                    trial['executions'][variant]=result;save(root/'summary.json',summary)
                for variant,_,_ in variants:
                    trial['comparisons'][variant]=compare(directory/'cpu',directory/variant,1e-8,2e-6)
                    for key in ('cluster_listing_sha256','cluster_counts_by_dimension'):
                        if trial['executions'][variant][key]!=trial['executions']['cpu'][key]:raise RuntimeError(f'{name}/{variant} {key} changed')
                trial['historical_cpu']=compare(source.parent,directory/'cpu',1e-8,2e-6)
                print('PASS',name,'repeat',repeat,{k:round(v['wall_s'],3) for k,v in trial['executions'].items()},flush=True)
            record['median_s']={v:statistics.median(t['executions'][v]['wall_s'] for t in record['trials']) for v,_,_ in variants}
            save(root/'summary.json',summary)
        summary['status']='passed';save(root/'summary.json',summary)
    except BaseException as error:
        summary.update(status='failed',error=f'{type(error).__name__}: {error}');save(root/'summary.json',summary);raise


if __name__=='__main__':main()
