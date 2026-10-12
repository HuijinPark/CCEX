#!/usr/bin/env python3
"""Check backend selection, unavailable-device handling and strict input validation."""
import argparse
import json
import os
from pathlib import Path
import subprocess

from benchmark import REPO, INPUTS, compare


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpu',type=Path,required=True)
    parser.add_argument('--cuda',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output=args.output.resolve()
    args.output.mkdir(parents=True,exist_ok=False)
    case=INPUTS/'Diamond_NV_natab_gcce'
    original=json.loads((case/'inputfiles/ccein_nstate0.json').read_text())
    original.update(nstep=2,rbath=10,order=1,profile=True)
    original['gyrofile']=str((case/original['gyrofile']).resolve())
    original['bathfile']=[str((case/p).resolve()) for p in original['bathfile']]

    def run(name,binary,overrides,expected=None,hidden=False,extra=(),
            reference='baseline',actual_backend='cpu'):
        directory=args.output/name
        directory.mkdir()
        config=dict(original,**overrides)
        config['outfile']=str(directory/'coherence')
        path=directory/'input.json'
        path.write_text(json.dumps(config))
        env=dict(os.environ,I_MPI_FABRICS='shm')
        if hidden: env['CUDA_VISIBLE_DEVICES']=''
        with (directory/'run.log').open('w') as log:
            result=subprocess.run([str(binary.resolve()),'-f',str(path),*extra],cwd=REPO,
                                  env=env,stdout=log,stderr=subprocess.STDOUT,timeout=30)
        output=(directory/'run.log').read_text()
        if expected:
            assert result.returncode!=0 and expected in output, (name,result.returncode,output[-1000:])
        else:
            assert result.returncode==0, (name,output[-1000:])
            if reference and name!=reference:
                tolerance=(1e-8,2e-6) if actual_backend=='cuda' else (1e-11,1e-11)
                compare(args.output/reference,directory,*tolerance)
            assert 'backend='+actual_backend in output, name
        print('PASS '+name,flush=True)

    run('baseline',args.cpu,{'backend':'cpu'})
    run('cpu-auto',args.cpu,{'backend':'auto'})
    run('cpu-reject-cuda',args.cpu,{'backend':'cuda'},'no CUDA support')
    run('invalid-backend',args.cpu,{'backend':True},'backend must be')
    run('invalid-dimension',args.cpu,{'gpu_options':{'min_dim':0}},'positive 32-bit integer')
    run('fractional-batch',args.cpu,{'gpu_options':{'evolution_batch_steps':1.5}},'positive 32-bit integer')
    run('string-batch',args.cpu,{'gpu_options':{'evolution_batch_steps':'64'}},'finite JSON number')
    run('invalid-profile',args.cpu,{'profile':'true'},'must be a boolean')
    run('invalid-refinement',args.cpu,{'gcce_refinement':'false'},'must be a boolean')
    run('zero-refinement-sweeps',args.cpu,{'gcce_refinement_max_sweeps':0},'positive 32-bit integer')
    run('fractional-refinement-sweeps',args.cpu,{'gcce_refinement_max_sweeps':1.5},'positive 32-bit integer')
    run('negative-detail-limit',args.cpu,{'pcce_detail_limit':-1},'nonnegative 32-bit integer')
    run('extended-refinement',args.cpu,{'gcce_refinement_max_sweeps':64})
    run('cpu-refinement-off',args.cpu,{'gcce_refinement':False},reference=None)
    run('cuda-refinement-off',args.cuda,{'backend':'cuda','gcce_refinement':False},reference='cpu-refinement-off',actual_backend='cuda')
    run('negative-cache',args.cpu,{'gpu_options':{'cache_mb':-1}},'nonnegative 32-bit integer')
    run('disabled-cache',args.cpu,{'gpu_options':{'cache_mb':0}})
    run('invalid-cluster-batch',args.cpu,{'gpu_options':{'evolution_batch_size':0}},'positive 32-bit integer')
    run('invalid-cluster-budget',args.cpu,{'gpu_options':{'cluster_batch_mb':0}},'positive 32-bit integer')
    run('invalid-partition',args.cpu,{'pcce_partition':'unknown'},'pcce_partition must be')
    run('invalid-connections',args.cpu,{'pcce_connections':'unknown'},'pcce_connections must be')
    run('cuda-hidden-cpu',args.cuda,{'backend':'cpu'},hidden=True)
    run('cuda-hidden-auto',args.cuda,{'backend':'auto'},hidden=True)
    run('cuda-hidden-required',args.cuda,{'backend':'cuda'},'CUDA initialization failed',hidden=True)
    run('cuda-reject-expm',args.cuda,{'backend':'cuda','propagator':'expm'},'supports method=CCE')
    run('cuda-reject-cli-dsj',args.cuda,{'backend':'cuda'},'supports method=CCE',extra=('-m','dsj'))
    run('cce-cpu',args.cpu,{'backend':'cpu','method':'CCE'},reference=None)
    run('cce-cuda-cli',args.cuda,{'backend':'cuda'},extra=('-m','CCE'),reference='cce-cpu',actual_backend='cuda')
    run('cce-auto-cpu',args.cuda,{'backend':'auto','method':'CCE','gpu_options':{'min_dim':256}},reference='cce-cpu')


if __name__=='__main__':
    main()
