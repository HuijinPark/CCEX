#!/usr/bin/env python3
# HISTORICAL (2026-09-12..16): option names and binaries of that time; does not run with the current build.
# See README.md here; the current check is ../check_adaptive/verify.sh.
"""Validate adaptive groups on full NV/PL6 bath fixtures across MPI/backends."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
REPO=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(REPO/'test/gpu_verification/check_gpu'))
from check_pcce import run_case
from benchmark import compare


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fixtures',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ranks',type=int,default=4)
    p.add_argument('--gpu-ranks',type=int,default=2)
    a=p.parse_args();f=a.fixtures.resolve();r=a.output.resolve();r.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((f/'search-manifest.json').read_text())
    for name,want in manifest['files'].items():assert hashlib.sha256((f/name).read_bytes()).hexdigest()==want
    for item in manifest['realizations']:assert hashlib.sha256((f/item['name']).read_bytes()).hexdigest()==item['sha256']
    families={family['name']:family for family in manifest['families']}
    cases=[('nv13c-cce',families['nv13c'],'nv13c-bath','cce',222),
           ('nv13c-gcce',families['nv13c'],'nv13c-bath','gcce',222)]
    cases += [(f'pl6-bath{i}',families['pl6'],f'pl6-bath-realization-{i}','cce',n) for i,n in [(1,788),(2,765),(3,768)]]
    result=dict(status='running',note='Functional comparisons, not controlled speed benchmarks.',cases=[])
    def save():(r/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    try:
        for name,family,bath,method,n in cases:
            cfg=dict(family['config'])
            for key in list(cfg):
                if key.startswith('pcce_') or key in ['sK','max_trial','max_iter','kmeans_pp','iter_detail']:cfg.pop(key)
            for key in ['gyrofile','qubitfile']:
                if cfg.get(key):cfg[key]=str(f/cfg[key])
            cfg.update(bathfile=[str(f/bath)],method=method,cluster_mode='adaptive',order=2,npulse=4,
                       nstep=33,deltat=family['end_ms']/32,nstate=0,adaptive_group_max_dim=4,
                       adaptive_cluster_max_dim=48,adaptive_merge_cutoff=1.0,profile=True)
            case=dict(case=name,config=cfg,selected_sites=n,executions=[]);result['cases'].append(case);save()
            reference=None
            for backend,ranks in [('cpu',a.ranks),('cuda',a.gpu_ranks),('auto',a.ranks)]:
                directory=r/name/backend;binary=REPO/f"bin/main.{'cpu' if backend=='cpu' else 'cuda'}.out"
                run,log=run_case(dict(cfg,backend=backend),directory,binary,ranks,600)
                assert run['returncode']==0,(name,backend,log[-3000:])
                groups=[(int(d),list(map(int,sites.split(',')))) for d,sites in
                        re.findall(r'\[adaptive-group\] id=\d+ dim=(\d+) sites=([\d,]+)',log)]
                assert sorted(i for _,g in groups for i in g)==list(range(n)),name
                assert all(d<=cfg['adaptive_group_max_dim'] for d,_ in groups)
                assert all(int(d)<=cfg['adaptive_cluster_max_dim'] for d in run['cluster_counts_by_dimension'])
                record=dict(backend=backend,ranks=ranks,execution=run,groups=groups,binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest())
                if reference:
                    record['comparison']=compare(reference[0],directory,1e-8,2e-6)
                    assert reference[1]==groups and reference[2]==run['cluster_listing_sha256']
                else:reference=(directory,groups,run['cluster_listing_sha256'])
                case['executions'].append(record);save();print('PASS',name,backend,'groups',len(groups),flush=True)
        result['status']='passed';save()
    except BaseException as error:
        result.update(status='failed',error=str(error));save();raise


if __name__=='__main__':main()
