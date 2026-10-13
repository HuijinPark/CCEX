#!/usr/bin/env python3
"""Check backend selection, unavailable-device handling and strict input validation."""
import argparse
import json
import os
import re
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

    # hidden: True hides every GPU (CUDA_VISIBLE_DEVICES empty), a string sets CUDA_VISIBLE_DEVICES to it; actual_backend
    # None: no "backend=" profile line is asked for (the packed kernel prints none), the GPU tolerance is used
    def run(name,binary,overrides,expected=None,hidden=False,extra=(),
            reference='baseline',actual_backend='cpu',env_extra=None):
        directory=args.output/name
        directory.mkdir()
        config=dict(original,**overrides)
        config['outfile']=str(directory/'coherence')
        path=directory/'input.json'
        path.write_text(json.dumps(config))
        env=dict(os.environ,I_MPI_FABRICS='shm')
        if hidden: env['CUDA_VISIBLE_DEVICES']=hidden if isinstance(hidden,str) else ''
        env.update(env_extra or {})
        with (directory/'run.log').open('w') as log:
            result=subprocess.run([str(binary.resolve()),'-f',str(path),*extra],cwd=REPO,
                                  env=env,stdout=log,stderr=subprocess.STDOUT,timeout=30)
        output=(directory/'run.log').read_text()
        if expected:
            assert result.returncode!=0 and expected in output, (name,result.returncode,output[-1000:])
        else:
            assert result.returncode==0, (name,output[-1000:])
            if reference and name!=reference:
                tolerance=(1e-11,1e-11) if actual_backend=='cpu' else (1e-8,2e-6)
                compare(args.output/reference,directory,*tolerance)
            if actual_backend: assert 'backend='+actual_backend in output, name
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
    # packed_max_dim 0 here: with the default (4) the GPU may take every cluster once it is up, and then there is no CPU
    # profile line and the tolerance is the GPU's (review 20261012-233058-452278 P3-05); the default is cce-auto-packed
    run('cce-auto-cpu',args.cuda,{'backend':'auto','method':'CCE','gpu_options':{'min_dim':256,'packed_max_dim':0}},
        reference='cce-cpu')
    run('cce-auto-packed',args.cuda,{'backend':'auto','method':'CCE','gpu_options':{'min_dim':256}},reference='cce-cpu',
        actual_backend=None)
    output=(args.output/'cce-auto-packed'/'run.log').read_text()
    assert '[gpu-packed] on' in output and re.search(r'\[gpu-packed\] clusters on the GPU \d+, on the CPU while the GPU started \d+',
                                                     output), 'cce-auto-packed'
    print('PASS cce-auto-packed log: the packed kernel on, its counts printed',flush=True)
    # P3-01: a card that cannot be opened - auto leaves every cluster on the CPU (the packed kernel off), cuda refuses, and
    # an explicit packed_max_dim with backend cpu (its only GPU work) refuses too
    for value in ('-1','99'):
        tag=value.replace('-','m')
        run('cce-hidden-auto-'+tag,args.cuda,{'backend':'auto','method':'CCE'},hidden=value,reference='cce-cpu')
        assert '[gpu-packed] off: the GPU could not be opened' in (args.output/('cce-hidden-auto-'+tag)/'run.log').read_text(), tag
        print('PASS cce-hidden-auto-'+tag+' log: the packed kernel off, every cluster on the CPU',flush=True)
        run('cce-hidden-cuda-'+tag,args.cuda,{'backend':'cuda','method':'CCE'},'CUDA initialization failed',hidden=value)
        run('cce-hidden-cpu-packed-'+tag,args.cuda,{'backend':'cpu','method':'CCE','gpu_options':{'packed_max_dim':4}},
            'Error: the packed GPU kernel',hidden=value)
    # P1-05: with the packed kernel on, clusters of a dimension it does not take (8 here) keep the ordinary CUDA batching
    # (warp_max_dim 0 here and in the out-of-memory cases below: the warp kernel, on by default since A40, would take them)
    third={'method':'CCE','order':3,'rbath':15}
    run('cce3-cpu',args.cpu,dict(third,backend='cpu'),reference=None)
    run('cce3-auto-packed-batch',args.cuda,dict(third,backend='auto',gpu_options={'min_dim':8,'evolution_batch_size':8,'warp_max_dim':0}),
        reference='cce3-cpu',actual_backend=None)
    output=(args.output/'cce3-auto-packed-batch'/'run.log').read_text()
    assert '[gpu-packed] on' in output and re.search(r'\[gpu-batch\] rank=\d+ dim=8 batches=[1-9]',output), 'cce3-auto-packed-batch'
    print('PASS cce3-auto-packed-batch log: dimension 8 in the ordinary CUDA batches beside the packed kernel',flush=True)
    # out of GPU memory (user 2026-10-13, A37 job 699): CCEX_TEST_GPU_OOM_BYTES=512 fails every device allocation of
    # 512 bytes or more, as a full card would. auto: those clusters on the CPU (CPU tolerance against the CPU build), the
    # run goes on; cuda: refused. cCCE and gCCE, in cluster batches (evolution_batch_size 8, the default) and one by one.
    oom={'CCEX_TEST_GPU_OOM_BYTES':'512'}
    def fell_back(name):
        output=(args.output/name/'run.log').read_text()
        assert '[gpu-memory-fallback] rank=0 dim=' in output and 'backend=cuda' not in output, name
        print('PASS '+name+' log: out of GPU memory, every cluster on the CPU',flush=True)
    for size in (8,1):
        name=f'cce3-oom-auto-batch{size}'
        run(name,args.cuda,dict(third,backend='auto',gpu_options={'min_dim':8,'packed_max_dim':0,'warp_max_dim':0,'evolution_batch_size':size}),
            reference='cce3-cpu',env_extra=oom)
        fell_back(name)
        name=f'gcce-oom-auto-batch{size}'
        run(name,args.cuda,{'backend':'auto','gpu_options':{'min_dim':6,'warp_max_dim':0,'evolution_batch_size':size}},env_extra=oom)
        fell_back(name)
    run('cce3-oom-cuda',args.cuda,dict(third,backend='cuda',gpu_options={'packed_max_dim':0,'warp_max_dim':0}),'out of memory',env_extra=oom)
    run('gcce-oom-cuda',args.cuda,{'backend':'cuda','gpu_options':{'warp_max_dim':0}},'out of memory',env_extra=oom)

    # evolution_update per path in the GPU build (A36, 2026-10-12): without the key the clusters on the CPU and the packed
    # kernel run recursive and the GPU large-cluster path exact, and the log says which; a given recursive is refused
    # where a cluster may reach the GPU large-cluster path (backend cuda, or auto), as before A36
    def chosen(name,value,why):
        output=(args.output/name/'run.log').read_text()
        assert re.search(r'^ +evolution_update +: +'+value+r'\b',output,re.M), (name,value)
        assert re.search(r'^ +chosen by +: +default \('+re.escape(why)+r'\)',output,re.M), (name,why)
        print('PASS '+name+' log: '+value+' ('+why+')',flush=True)
    chosen('cuda-refinement-off','exact','backend=cuda: every cluster on the GPU, exact there')
    chosen('cuda-hidden-auto','recursive','clusters on the CPU: recursive; clusters on the GPU (backend=auto, warp kernel included): exact')
    chosen('cce-cuda-cli','recursive','packed kernel: recursive; clusters on the GPU (backend=cuda, warp kernel included): exact')
    chosen('cce-auto-cpu','recursive','clusters on the CPU: recursive; clusters on the GPU (backend=auto, warp kernel included): exact')
    chosen('cce-auto-packed','recursive','clusters on the CPU and packed kernel: recursive; clusters on the GPU (backend=auto, warp kernel included): exact')
    chosen('baseline','recursive','supported path')
    run('gpu-build-cpu-recursive',args.cuda,{'backend':'cpu','evolution_update':'recursive'})
    run('auto-reject-recursive',args.cuda,{'backend':'auto','evolution_update':'recursive'},'this run has backend=auto')
    run('cuda-reject-recursive',args.cuda,{'backend':'cuda','evolution_update':'recursive'},'this run has backend=cuda')
    run('cce-cuda-reject-recursive',args.cuda,{'backend':'cuda','evolution_update':'recursive'},'this run has backend=cuda',
        extra=('-m','CCE'))

    # the warp kernel (A40, gpu_options.warp_max_dim; user 2026-10-13: default 27 in the GPU build with backend auto or
    # cuda and method CCE or gCCE, 0 otherwise): a given value is refused where it cannot run; the log says the value and why
    def warp(name,value,why):
        output=(args.output/name/'run.log').read_text()
        assert re.search(r'^ +gpu_options.warp_max_dim *: +'+str(value)+r'\b.*\n +chosen by +: +'+re.escape(why),output,re.M), \
            (name,value,why)
        print('PASS '+name+' log: warp_max_dim '+str(value)+' ('+why+')',flush=True)
    warp('baseline',0,'default (CPU build)')
    warp('cuda-hidden-cpu',0,'default (backend=cpu)')
    warp('cuda-hidden-auto',27,'default (backend=auto, method=gCCE)')
    warp('cuda-refinement-off',27,'default (backend=cuda, method=gCCE)')
    warp('cce-cuda-cli',27,'default (backend=cuda, method=CCE)')
    run('warp-reject-cpu-build',args.cpu,{'backend':'auto','gpu_options':{'warp_max_dim':27}},'needs the GPU binary')
    run('warp-reject-backend-cpu',args.cuda,{'backend':'cpu','gpu_options':{'warp_max_dim':27}},'works with backend auto or cuda')
    run('warp-reject-dimension',args.cuda,{'backend':'auto','gpu_options':{'warp_max_dim':5}},'has no dimension up to it')
    run('warp-reject-form',args.cuda,{'backend':'auto','gpu_options':{'warp_max_dim':33}},'an integer up to 32')
    run('warp-reject-old-key',args.cuda,{'backend':'auto','gpu_options':{'warp_dims':[6]}},'unknown key gpu_options.warp_dims')
    run('warp-reject-pcce',args.cuda,{'backend':'auto','gpu_options':{'warp_max_dim':27}},'works with method=CCE or gCCE',
        extra=('-m','pCCE'))
    # (this input has no sK, so the pCCE run stops at clustering - after the log has printed the settings)
    run('warp-default-pcce',args.cuda,{'backend':'auto'},'pCCE sK',extra=('-m','pCCE'))
    warp('warp-default-pcce',0,'default (not method=CCE (coherence) or gCCE)')
    run('warp-off-gcce',args.cuda,{'backend':'auto','gpu_options':{'warp_max_dim':0}},reference=None)
    warp('warp-off-gcce',0,'given')
    # the default switches itself off where the kernel does not apply (here gcce_refinement), and the run goes on
    run('warp-default-refinement',args.cuda,{'backend':'auto','gcce_refinement':True},reference=None)
    assert '[gpu-warp] off: gcce_refinement' in (args.output/'warp-default-refinement'/'run.log').read_text()
    print('PASS warp-default-refinement log: [gpu-warp] off',flush=True)

if __name__=='__main__':
    main()
