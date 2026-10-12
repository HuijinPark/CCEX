#!/usr/bin/env python3
# HISTORICAL (2026-09-12..16): option names and binaries of that time; does not run with the current build.
# See README.md here; the current check is ../check_adaptive/verify.sh.
"""Check adaptive CCE/GCCE against existing full-bath and subset calculations."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import re
import sys

REPO=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(REPO/'test/gpu_verification/check_gpu'))
from check_pcce import run_case
from benchmark import compare,values


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cpu',type=Path,default=REPO/'bin/main.cpu.out')
    p.add_argument('--cuda',type=Path,default=REPO/'bin/main.cuda.out')
    p.add_argument('--cpu-only',action='store_true')
    a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    base=json.loads((REPO/'test/cluster_verification/check_bundle/ensemble.json').read_text())
    summary=dict(status='running',records=[],oracles=[],accuracy=[])
    def save(): (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    def run(name,cfg,ranks=1,backend='cpu',old=False,reference=None,error=None,extra=()):
        binary=(a.baseline if old else a.cpu if backend=='cpu' else a.cuda).resolve()
        directory=root/name
        execution,log=run_case(dict(cfg,backend=backend),directory,binary,ranks,180,extra=extra)
        record=dict(case=name,backend=backend,ranks=ranks,baseline=old,execution=execution,
                    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest())
        if error:
            assert execution['returncode']!=0 and error in log,(name,log[-3000:]);record['expected_error']=error
        else:
            assert execution['returncode']==0,(name,log[-4000:])
            if cfg.get('cluster_mode')=='adaptive':
                groups=[]
                for dim,ids in re.findall(r'\[adaptive-group\] id=\d+ dim=(\d+) sites=([\d,]+)',log):
                    groups.append(dict(dim=int(dim),sites=list(map(int,ids.split(',')))))
                sites=[i for g in groups for i in g['sites']]
                assert sorted(sites)==list(range(len(sites))),'omitted or repeated site'
                record['groups']=groups
                assert all(int(d)<=cfg['adaptive_cluster_max_dim'] for d in execution['cluster_counts_by_dimension'])
            if reference:
                record['comparison']=compare(reference,directory,1e-8,2e-6)
        summary['records'].append(record);save();print('PASS',name,flush=True)
        return directory,record
    def standard(cfg,order):
        return {**{k:v for k,v in cfg.items() if not k.startswith('adaptive_') and k!='cluster_mode'},'order':order}
    def curve(directory):
        files=list(directory.glob('*_wiDiv'));assert len(files)==1
        return [row[1] for row in values(files[0])]
    try:
        # One or two groups covering the full bath: the only propagated term is
        # the full Hamiltonian, independently selected by the old full-order CCE.
        for method,state,pulse in itertools.product(['cce','gcce'],[0,1],[1,4]):
            name=f'full-{method}-s{state}-p{pulse}'
            cfg=dict(base,method=method,nstate=state,npulse=pulse,adaptive_group_max_dim=16,
                     adaptive_cluster_max_dim=64 if method=='cce' else 192,adaptive_merge_cutoff=0)
            reference,_=run(name+'-old',standard(cfg,6),old=True)
            _,record=run(name+'-cpu',cfg,reference=reference)
            assert len(record['groups'])==2
            if not a.cpu_only:run(name+'-cuda',cfg,2,'cuda',reference=reference)
            if pulse==1:
                run(name+'-one-group',dict(cfg,order=1,adaptive_group_max_dim=64),4,reference=reference)
        # Singletons reproduce conventional CCE with the same pair graph,
        # including random-state generation and gCCE's nontrivial L0 factor.
        for method,state in itertools.product(['cce','gcce'],[0,1]):
            cfg=dict(base,method=method,nstate=state,npulse=4,adaptive_group_max_dim=2,
                     adaptive_cluster_max_dim=12,adaptive_merge_cutoff=0)
            name=f'singletons-{method}-s{state}'
            ref,_=run(name+'-old',standard(cfg,2),old=True)
            _,rec=run(name+'-mpi4',cfg,4,reference=ref);assert len(rec['groups'])==6
            if not a.cpu_only:run(name+'-auto',dict(cfg,gpu_min_dim=8),2,'auto',reference=ref)
        # Independent subset oracle for unequal 3/2/1 groups. Group triplets and
        # pair unions occupy the same site-count bucket; both must be retained.
        bathrows=(REPO/'test/cluster_verification/check_bundle/bath_6').read_text().splitlines()[1:]
        for method,pulse in itertools.product(['cce','gcce'],[1,4]):
            name=f'unequal-{method}-p{pulse}'
            cfg=dict(base,method=method,npulse=pulse,adaptive_cluster_max_dim=192)
            ref,rec=run(name+'-cpu',cfg)
            assert [len(g['sites']) for g in rec['groups']]==[3,2,1]
            _,mpi=run(name+'-mpi4',cfg,4,reference=ref);assert rec['groups']==mpi['groups']
            if not a.cpu_only:
                _,gpu=run(name+'-cuda',cfg,2,'cuda',reference=ref);assert rec['groups']==gpu['groups']
            groups=[g['sites'] for g in rec['groups']];signals={}
            zero=[1+0j]*cfg['nstep']
            if method=='gcce':
                z,_=run(name+'-zero-old',standard(cfg,0),old=True);zero=curve(z)
            expected=list(zero)
            for order in [1,2]:
                for chosen in itertools.combinations(range(3),order):
                    ids=sorted(i for g in chosen for i in groups[g])
                    f=root/(name+'-sites-'+'-'.join(map(str,chosen)))
                    f.write_text(f'{len(ids)+1}.0 0 0 0\n'+'\n'.join(bathrows[i] for i in ids)+'\n')
                    subset=dict(standard(cfg,len(ids)),bathfile=[str(f)])
                    d,_=run(name+'-subset-'+'-'.join(map(str,chosen)),subset,old=True)
                    signal=curve(d);signals[chosen]=signal
                    for t in range(len(zero)):
                        if order==1:expected[t]*=signal[t]/zero[t]
                        else:expected[t]*=signal[t]*zero[t]/(signals[(chosen[0],)][t]*signals[(chosen[1],)][t])
            got=curve(ref);error=max(abs(x-y) for x,y in zip(got,expected))
            assert error<2e-8,(name,error)
            summary['oracles'].append(dict(case=name,independent_subset_product_max_abs_error=error));save()
            exact,_=run(name+'-exact-old',standard(cfg,6),old=True)
            conventional,_=run(name+'-cce2-old',standard(cfg,2),old=True)
            exactcurve=curve(exact)
            summary['accuracy'].append(dict(case=name,
                adaptive_max_abs_error=max(abs(x-y) for x,y in zip(got,exactcurve)),
                conventional_order2_max_abs_error=max(abs(x-y) for x,y in zip(curve(conventional),exactcurve)),
                note='Sampled output error against full six-spin Hamiltonian, not an error bound for larger baths.'))
            save()
        # Mixed spin quantum numbers: 2,3,6 local dimensions, pair dimensions 6,12,18.
        gyro=root/'mixed.gyro';gyro.write_text('A 0.5 6.728284\nB 1 -2.71261804\nC 2.5 -1.784\n')
        bath=root/'mixed.bath';bath.write_text('4.0 0 0 0\n5 1 2 A\n6 2 1 B\n1 7 3 C\n')
        for method in ['cce','gcce']:
            cfg=dict(base,method=method,gyrofile=str(gyro),bathfile=[str(bath)],adaptive_group_max_dim=6,
                     adaptive_cluster_max_dim=54,adaptive_merge_cutoff=1e100)
            ref,_=run('mixed-'+method+'-old',standard(cfg,2),old=True)
            _,rec=run('mixed-'+method+'-cpu',cfg,reference=ref)
            assert [g['dim'] for g in rec['groups']]==[2,3,6]
            if not a.cpu_only:run('mixed-'+method+'-cuda',cfg,2,'cuda',reference=ref)
        # P1 sites each carry an electron AND an added 15N spin (site dimension 4).
        p1=json.loads((REPO/'test/gpu_verification/check_gpu/examples/p1_gcce4_single.json').read_text())
        source=REPO/p1['bathfile'][0];rows=source.read_text().splitlines()[1:3]
        bath=root/'p1-two.bath';bath.write_text('3.0 0 0 0\n'+'\n'.join(rows)+'\n')
        p1.update(bathfile=[str(bath)],cluster_mode='adaptive',adaptive_group_max_dim=16,
                  adaptive_cluster_max_dim=48,adaptive_merge_cutoff=0,order=1,nstep=9,npulse=4)
        ref,_=run('p1-old',standard(p1,2),old=True)
        _,rec=run('p1-adaptive-cpu',p1,reference=ref);assert [g['dim'] for g in rec['groups']]==[16]
        if not a.cpu_only:run('p1-adaptive-cuda',p1,2,'cuda',reference=ref)
        run('p1-site-cap',dict(p1,adaptive_group_max_dim=2),error='indivisible site exceeds')
        # Existing pCCE behavior and its cache identity are independent of this mode.
        oldpcce=json.loads((REPO/'test/gpu_verification/check_gpu/examples/pcce_ensemble.json').read_text())
        ref,_=run('pcce-old',oldpcce,old=True)
        run('pcce-unchanged',oldpcce,4,reference=ref)
        if not a.cpu_only:run('pcce-unchanged-cuda',oldpcce,2,'cuda',reference=ref)
        bad=[('order0',dict(order=0),'positive integer group order'),
             ('fractional-cap',dict(adaptive_group_max_dim=3.5),'positive 32-bit'),
             ('boolean-cap',dict(adaptive_group_max_dim=True),'must'),
             ('negative-cutoff',dict(adaptive_merge_cutoff=-1),'nonnegative'),
             ('site-cap',dict(adaptive_group_max_dim=1),'indivisible site exceeds'),
             ('pair-cap',dict(adaptive_cluster_max_dim=3),'required pair'),
             ('qubit-cap',dict(method='gcce',adaptive_cluster_max_dim=1),'exceeds adaptive_cluster_max_dim'),
             ('pcce-method',dict(method='pcce'),'pCCE requires sK'),
             ('fixed-size-option',dict(sK=2),'not an adaptive cluster option'),
             ('cache-option',dict(pcce_cache_mode='read'),'does not use pcce_cache_mode'),
             ('invalid-mode',dict(cluster_mode='magic'),'standard or adaptive'),
             ('standard-with-options',dict(cluster_mode='standard'),'require cluster_mode=adaptive'),
             ('limited-nk',dict(nk=[0,0,1]),'requires nk=0'),
             ('distance-cut',dict(rdipcut=base['rdip']),'rdipcut < rdip')]
        for name,change,error in bad:run('reject-'+name,dict(base,**change),2,error=error)
        cfg=dict(base);del cfg['adaptive_cluster_max_dim']
        run('missing-adaptive_cluster_max_dim',cfg,2,error='adaptive requires adaptive_cluster_max_dim')
        for key,default in [('adaptive_group_max_dim',base['adaptive_cluster_max_dim']),
                            ('adaptive_merge_cutoff',0)]:
            ref,_=run('explicit-default-'+key,dict(base,**{key:default}))
            cfg=dict(base);del cfg[key]
            run('missing-'+key,cfg,2,reference=ref)
        run('cli-method-override',base,2,error='requires method=cce or gcce',extra=('-m','pcce'))
        summary['status']='passed';save()
    except BaseException as error:
        summary.update(status='failed',error=f'{type(error).__name__}: {error}');save();raise


if __name__=='__main__':main()
