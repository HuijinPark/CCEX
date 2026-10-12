#!/usr/bin/env python3
# HISTORICAL (2026-09-12..16): option names and binaries of that time; does not run with the current build.
# See README.md here; the current check is ../check_adaptive/verify.sh.
"""End-to-end refinement checks against complete small-bath calculations."""
import argparse
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys

REPO=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(REPO/'test/gpu_verification/check_gpu'))
from check_pcce import run_case
from benchmark import values


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cpu',type=Path,default=REPO/'bin/main.cpu.out')
    p.add_argument('--cuda',type=Path)
    p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    base=json.loads((REPO/'test/cluster_verification/check_bundle/ensemble.json').read_text())
    report=dict(status='running',cases=[])
    def save(): (root/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    def run(name,cfg,expected='validated_full_bath_model',rounds=4,ranks=1,backend='cpu',oracle=True,extra=()):
        inputfile=root/(name+'.json');inputfile.write_text(json.dumps(dict(cfg,backend=backend)))
        binary=(a.cuda if backend=='cuda' else a.cpu).resolve();dest=root/name
        command=[sys.executable,str(REPO/'test/cluster_verification/check_adaptive/oracle/adaptive_refine.py'),'--input',str(inputfile),
                 '--binary',str(binary),'--output',str(dest),'--atol','1e-10','--rtol','1e-9',
                 '--max-rounds',str(rounds),'--ranks',str(ranks),*extra]
        with (root/(name+'.driver.log')).open('w') as log:
            proc=subprocess.run(command,cwd=REPO,stdout=log,stderr=subprocess.STDOUT,timeout=900)
        result=json.loads((dest/'summary.json').read_text())
        assert result['status']==expected,(name,result.get('status'),result.get('error'))
        assert proc.returncode==(0 if expected=='validated_full_bath_model' else 2)
        rec=dict(name=name,backend=backend,ranks=ranks,status=result['status'],rounds=len(result['rounds']),
                 final_groups=result['final_groups'],directory=str(dest))
        ids=sorted(s for g in rec['final_groups'] for s in g)
        assert ids==list(range(len(ids)))
        oracle_errors=[row['whole_bath_error']['scaled_error'] for row in result['rounds'] if 'whole_bath_error' in row]
        assert all(y<=x+1e-7 for x,y in zip(oracle_errors,oracle_errors[1:])),(name,oracle_errors)
        if oracle:
            full={k:v for k,v in cfg.items() if k!='cluster_mode' and not k.startswith('adaptive_')}
            full.update(order=len(ids),rdip=1e6,backend='cpu')
            reference=root/(name+'-complete')
            execution,log=run_case(full,reference,a.baseline.resolve(),ranks=1,timeout=300)
            assert execution['returncode']==0,(name,log[-2000:])
            differences=[]
            for original in reference.glob('*_wiDiv*'):
                candidate=Path(result['final_directory'])/original.name
                before=values(original);after=values(candidate);assert len(before)==len(after)
                for (_,x),(_,y) in zip(before,after):
                    assert math.isfinite(abs(y)) and abs(x-y)<=1e-8+2e-6*abs(x),(name,original.name,x,y)
                    differences.append(abs(x-y))
            assert differences,name
            rec['complete_comparison']=dict(max_abs_error=max(differences),points=len(differences),
                observable='wiDiv coherence; noDiv depends on partition and is not the complete-bath oracle')
        report['cases'].append(rec);save();print('PASS',name,rec['status'],rec['final_groups'],flush=True)
        return result
    for method,state,pulses in itertools.product(('cce','gcce'),(0,2),(1,4,16)):
        cfg=dict(base,method=method,nstate=state,npulse=pulses,adaptive_cluster_max_dim=192 if method=='gcce' else 64)
        name=f'{method}-s{state}-p{pulses}'
        result=run(name+'-cpu',cfg,ranks=2 if state else 1)
        if a.cuda:
            gpu=run(name+'-cuda',cfg,ranks=2,backend='cuda')
            assert result['final_groups']==gpu['final_groups'],name
    run('round-limit',base,expected='round_limit',rounds=0,oracle=False)
    run('probe-budget',base,expected='probe_budget_limited',oracle=False,extra=('--max-probes','1'))
    # Auto connection audit must also accept the explicit refined partition.
    auto=dict(base,adaptive_connections='auto',adaptive_connection_budget=0)
    auto.pop('rdip');run('auto-connections',auto)
    # Real dense h-BN coordinates and isotope identities; point-dipole fixture,
    # without the material's ab-initio hyperfine/quadrupole tensors.
    stored=REPO/'test/code_verification/bath_data/ensemble_data'
    q=list(map(float,(stored/'bath_h_BN_after_relax_Flat_1L_Defect').read_text().split()))
    rows=[line.split() for line in (stored/'bath_h_BN_after_relax_Flat_1L_1').read_text().splitlines()[1:]]
    nearest=sorted(enumerate(rows),key=lambda item:math.dist(q,list(map(float,item[1][:3]))))[:4]
    bath=root/'hbn-four.bath';bath.write_text('5 0 0 0\n'+'\n'.join(' '.join(row) for _,row in nearest)+'\n')
    report['hbn_fixture']=dict(original_sites=[i for i,_ in nearest],species=[r[3] for _,r in nearest],
        note='Four nearest stored nuclei; functional/convergence test, not a quantitative h-BN material prediction.')
    for method in ('cce','gcce'):
        cfg=dict(base,method=method,bathfile=[str(bath)],gyrofile=str(stored/'Gyro_h_BN'),
                 Qubit=dict(nqubit=1,qubit=[dict(name='VB',spin=1,gyro=-17608.59705,xyz=q,alphams=1,betams=0)]),
                 adaptive_group_max_dim=36,adaptive_cluster_max_dim=324 if method=='gcce' else 108,
                 adaptive_merge_cutoff=1e100,deltat=.0001,npulse=4)
        run('hbn-'+method,cfg,ranks=2)
        if a.cuda:run('hbn-'+method+'-cuda',cfg,ranks=2,backend='cuda')
    report['status']='passed';save()


if __name__=='__main__':main()
