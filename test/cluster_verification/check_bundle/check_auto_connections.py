#!/usr/bin/env python3
# HISTORICAL (2026-09-12..16): option names and binaries of that time; does not run with the current build.
# See README.md here; the current check is ../check_adaptive/verify.sh.
"""Auto-link integration: old full graph, independent subsets, MPI/CUDA, budgets."""
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
    p.add_argument('--cpu-only',action='store_true')
    p.add_argument('--fixtures',type=Path)
    p.add_argument('--real-family',choices=['nv13c','pl6'],help='Run only this fixture family; requires --fixtures')
    p.add_argument('--ranks',type=int,default=4)
    p.add_argument('--gpu-ranks',type=int,default=2)
    a=p.parse_args()
    if a.real_family and not a.fixtures:p.error('--real-family requires --fixtures')
    root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    base=json.loads((REPO/'test/cluster_verification/check_bundle/auto_ensemble.json').read_text())
    summary=dict(status='running',records=[],oracles=[],accuracy=[],real_baths=[])
    def save(): (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    def run(name,cfg,ranks=1,backend='cpu',old=False,reference=None,error=None):
        binary=a.baseline.resolve() if old else REPO/f"bin/main.{'cpu' if backend=='cpu' else 'cuda'}.out"
        directory=root/name
        execution,log=run_case(dict(cfg,backend=backend),directory,binary,ranks,600)
        rec=dict(case=name,backend=backend,ranks=ranks,baseline=old,execution=execution,
                 binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest())
        if error:
            assert execution['returncode']!=0 and error in log,(name,log[-4000:]);rec['expected_error']=error
        else:
            assert execution['returncode']==0,(name,log[-4000:])
            compare(directory,directory,1e-8,2e-6)
            if cfg.get('adaptive_connections')=='auto':
                m=re.search(r'\[adaptive-connections\] (.*)',log);assert m,name
                stats=dict(x.split('=') for x in m[1].split());rec['connections']=stats
                assert float(stats['omitted_score'])<=float(stats['budget'])+1e-14
                assert int(stats['positive_pairs'])==int(stats['kept_pairs'])+int(stats['omitted_pairs'])
                rec['groups']=[dict(dim=int(d),sites=list(map(int,ids.split(','))))
                    for d,ids in re.findall(r'\[adaptive-group\] id=\d+ dim=(\d+) sites=([\d,]+)',log)]
                rec['links']=[tuple(map(int,ids.split(','))) for ids in re.findall(r'\[adaptive-link\] groups=([\d,]+)',log)]
                assert len(rec['links'])==int(stats['kept_pairs'])
                ids=[i for g in rec['groups'] for i in g['sites']];assert sorted(ids)==list(range(len(ids)))
                assert all(int(d)<=cfg['adaptive_cluster_max_dim'] for d in execution['cluster_counts_by_dimension'])
            if reference:rec['comparison']=compare(reference,directory,1e-8,2e-6)
        summary['records'].append(rec);save();print('PASS',name,flush=True)
        return directory,rec
    def distance(cfg):
        return {**{k:v for k,v in cfg.items() if k not in ['adaptive_connections','adaptive_connection_budget']},'rdip':1e6}
    def standard(cfg,order):
        return {**{k:v for k,v in cfg.items() if not k.startswith('adaptive_') and k!='cluster_mode'},'order':order,'rdip':1e6}
    def curve(path):
        fs=list(path.glob('*_wiDiv'));assert len(fs)==1
        return [row[1] for row in values(fs[0])]
    def same_graph(x,y):
        assert x['groups']==y['groups'] and x['links']==y['links']
        assert x['connections']==y['connections']
    try:
        if not a.real_family:
            # Each backend uses the same graph for a whole curve; pulse/state changes
            # affect dynamics but not this coupling-norm selection at fixed tmax.
            fixed_graph=None
            for method,state,pulse in itertools.product(['cce','gcce'],[0,1],[1,4,16]):
                cfg=dict(base,method=method,nstate=state,npulse=pulse,
                         adaptive_cluster_max_dim=64 if method=='cce' else 192)
                name=f'{method}-s{state}-p{pulse}'
                old,_=run(name+'-old-full-graph',distance(cfg),old=True)
                run(name+'-zero-budget',dict(cfg,adaptive_connection_budget=0),reference=old)
                ref,rec=run(name+'-pruned-cpu',cfg)
                assert [len(g['sites']) for g in rec['groups']]==[3,2,1]
                assert len(rec['links'])==1 and int(rec['connections']['omitted_pairs'])==2
                if fixed_graph is None:fixed_graph=rec
                else:same_graph(fixed_graph,rec)
                if not a.cpu_only:
                    _,gpu=run(name+'-cuda',cfg,a.gpu_ranks,'cuda',reference=ref);same_graph(rec,gpu)
                if state==0 and pulse==4:
                    _,mpi=run(name+'-mpi',cfg,a.ranks,reference=ref);same_graph(rec,mpi)
                    if not a.cpu_only:
                        _,auto=run(name+'-auto',cfg,2,'auto',reference=ref);same_graph(rec,auto)
                    _,wide=run(name+'-wide-rdip',dict(cfg,rdip=1e6),reference=ref)
                    assert wide['groups']==rec['groups'] and wide['links']==rec['links']
                if state==0 and pulse in [1,16]:
                    # Independently assemble the selected graph's pair corrections
                    # from old full-Hamiltonian calculations for each subset.
                    rows=(REPO/'test/cluster_verification/check_bundle/bath_6').read_text().splitlines()[1:]
                    zero=[1+0j]*cfg['nstep']
                    if method=='gcce':z,_=run(name+'-zero-old',standard(cfg,0),old=True);zero=curve(z)
                    signals={};expected=list(zero)
                    subsets=[(g,) for g in range(len(rec['groups']))]+rec['links']
                    for chosen in subsets:
                        ids=sorted(i for g in chosen for i in rec['groups'][g]['sites'])
                        bath=root/(name+'-subset-'+'-'.join(map(str,chosen))+'.bath')
                        bath.write_text(f'{len(ids)+1}.0 0 0 0\n'+'\n'.join(rows[i] for i in ids)+'\n')
                        d,_=run(name+'-subset-'+'-'.join(map(str,chosen)),dict(standard(cfg,len(ids)),bathfile=[str(bath)]),old=True)
                        sig=curve(d);signals[tuple(chosen)]=sig
                        for t in range(len(zero)):
                            if len(chosen)==1:expected[t]*=sig[t]/zero[t]
                            else:expected[t]*=sig[t]*zero[t]/(signals[(chosen[0],)][t]*signals[(chosen[1],)][t])
                    got=curve(ref);err=max(abs(x-y) for x,y in zip(got,expected));assert err<2e-8
                    summary['oracles'].append(dict(case=name,max_abs_error=err))
                    exact,_=run(name+'-exact',standard(cfg,6),old=True)
                    summary['accuracy'].append(dict(case=name,max_abs_error=max(abs(x-y) for x,y in zip(got,curve(exact))),
                        note='Measured small-bath error, not guaranteed by the connection budget.'));save()
            for budget in [.001,.005,.05]:
                _,rec=run('budget-'+str(budget),dict(base,adaptive_connection_budget=budget))
                assert int(rec['connections']['kept_pairs'])=={.001:3,.005:2,.05:0}[budget]
            _,rec=run('longer-time',dict(base,deltat=.02))
            assert len(rec['links'])==2
            run('time-zero',dict(base,nstep=1,adaptive_connection_budget=0))
            run('one-group-order1',dict(base,order=1,adaptive_group_max_dim=64,adaptive_merge_cutoff=0,adaptive_connection_budget=0))
            _,cut=run('loose-radius-cut',dict(base,rdip=.1,adaptive_connection_budget=1e6))
            assert len(cut['groups'])==6 and not cut['links'] and int(cut['connections']['radius_excluded_pairs'])==15
            # Mixed I and defect subspins also participate in the all-distance graph.
            gyro=root/'mixed.gyro';gyro.write_text('A 0.5 6.728284\nB 1 -2.71261804\nC 2.5 -1.784\n')
            bath=root/'mixed.bath';bath.write_text('4.0 0 0 0\n5 1 2 A\n6 2 1 B\n1 7 3 C\n')
            for method in ['cce','gcce']:
                cfg=dict(base,method=method,gyrofile=str(gyro),bathfile=[str(bath)],adaptive_group_max_dim=6,
                    adaptive_cluster_max_dim=54,adaptive_merge_cutoff=1e100,adaptive_connection_budget=0)
                old,_=run('mixed-'+method+'-old',distance(cfg),old=True)
                ref,rec=run('mixed-'+method+'-cpu',cfg,reference=old)
                assert [g['dim'] for g in rec['groups']]==[2,3,6]
                if not a.cpu_only:run('mixed-'+method+'-cuda',cfg,2,'cuda',reference=ref)
            p1=json.loads((REPO/'test/gpu_verification/check_gpu/examples/p1_gcce4_single.json').read_text())
            rows=(REPO/p1['bathfile'][0]).read_text().splitlines()[1:3]
            bath=root/'p1-two.bath';bath.write_text('3.0 0 0 0\n'+'\n'.join(rows)+'\n')
            p1.update(bathfile=[str(bath)],cluster_mode='adaptive',adaptive_connections='auto',
                adaptive_group_max_dim=4,adaptive_cluster_max_dim=48,adaptive_merge_cutoff=0,
                adaptive_connection_budget=0,order=2,nstep=9,npulse=16);p1.pop('rdip')
            old,_=run('p1-old',distance(p1),old=True)
            ref,rec=run('p1-auto-cpu',p1,reference=old);assert [g['dim'] for g in rec['groups']]==[4,4]
            if not a.cpu_only:run('p1-auto-cuda',p1,2,'cuda',reference=ref)
            bad=[('narrow-rdip',dict(rdip=.1,adaptive_connection_budget=0),'budget exceeded'),
                 ('lower-distance-cut',dict(rdipcut=1000,adaptive_connection_budget=0),'budget exceeded'),
                 ('order1-small-budget',dict(order=1,adaptive_connection_budget=0),'budget exceeded'),
                 ('connection-mode',dict(adaptive_connections='magic'),'distance or auto'),
                 ('negative-budget',dict(adaptive_connection_budget=-1),'nonnegative'),
                 ('boolean-budget',dict(adaptive_connection_budget=True),'finite JSON number'),
                 ('zero-rdip',dict(rdip=0),'finite and positive'),
                 ('huge-rdip',dict(rdip=1e100),'finite and positive'),
                 ('string-rdip',dict(rdip='auto'),'finite JSON number'),
                 ('null-rdip',dict(rdip=None),'finite JSON number'),
                 ('budget-in-distance',dict(adaptive_connections='distance',rdip=30),'requires adaptive_connections=auto')]
            for name,change,error in bad:run('reject-'+name,dict(base,**change),2,error=error)
            cfg=dict(base);cfg.pop('adaptive_connection_budget');run('missing-budget',cfg,2,error='requires a finite nonnegative')
            cfg=distance(base);cfg.pop('rdip');run('distance-missing-rdip',cfg,2,error='rdip')
        if a.fixtures:
            fixture=a.fixtures.resolve();manifest=json.loads((fixture/'manifest.json').read_text())
            families={f['name']:f for f in manifest['families']}
            for family,method in [('nv13c','cce'),('nv13c','gcce'),('pl6','cce')]:
                if a.real_family and family!=a.real_family:continue
                cfg=dict(families[family]['config'])
                for key in ['qubitfile','gyrofile']:
                    if key in cfg:cfg[key]=str(fixture/cfg[key])
                cfg['bathfile']=[str(fixture/p) for p in cfg['bathfile']]
                for key in ['sK','max_trial','max_iter','kmeans_pp','iter_detail','pcce_partition','rdip']:cfg.pop(key,None)
                cfg.update(method=method,cluster_mode='adaptive',adaptive_connections='auto',
                    adaptive_group_max_dim=4,adaptive_cluster_max_dim=48,adaptive_merge_cutoff=1,
                    adaptive_connection_budget=.01,order=2,npulse=4,nstep=9,deltat=.0005,nstate=0)
                name='real-'+family+'-'+method
                # Short-time functional checks keep the old product's known
                # long-time PL6 singularities separate from graph validation.
                old,_=run(name+'-old-full-graph',distance(cfg),a.ranks,old=True)
                _,zero=run(name+'-zero',dict(cfg,adaptive_connection_budget=0),a.ranks,reference=old)
                ref,rec=run(name+'-pruned-cpu',cfg,a.ranks)
                if not a.cpu_only:
                    _,gpu=run(name+'-pruned-cuda',cfg,a.gpu_ranks,'cuda',reference=ref);same_graph(rec,gpu)
                summary['real_baths'].append(dict(case=name,config=cfg,groups=len(rec['groups']),
                    selected_sites=sum(len(g['sites']) for g in rec['groups']),full_pairs=len(zero['links']),
                    kept_pairs=len(rec['links']),connections=rec['connections'],
                    full_graph_difference=max(abs(x-y) for x,y in zip(curve(old),curve(ref)))));save()
        summary['status']='passed';save()
    except BaseException as error:
        summary.update(status='failed',error=f'{type(error).__name__}: {error}');save();raise

if __name__=='__main__':main()
