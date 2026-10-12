#!/usr/bin/env python3
"""Inventory actual clusters without invoking the dynamics routines."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'gpu_verification/check_gpu'))
from check_pcce import run_case
from check_baths import save,sha


def resolve(config,fixture):
    config=dict(config)
    for key in ['qubitfile','gyrofile']:
        if config.get(key):config[key]=str(fixture/config[key])
    config['bathfile']=[str(fixture/p) for p in config['bathfile']]
    return config


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fixtures',type=Path,required=True)
    p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ranks',type=int,default=4)
    a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    fixture=a.fixtures.resolve();manifest=json.loads((fixture/'manifest.json').read_text())
    for name,want in manifest['files'].items():assert sha(fixture/name)==want,name
    result=dict(status='running',ranks=a.ranks,cases=[])
    for family in manifest['families']:
        for size in [1,2,3]:
            groups=None
            for order in [1,2,3]:
                name=f"{family['name']}-o{order}-k{size}"
                cfg=dict(resolve(family['config'],fixture),sK=size,order=order,nstep=33,deltat=family['end_ms']/32,backend='cpu')
                execution,log=run_case(cfg,root/name,a.binary.resolve(),a.ranks,1800)
                assert execution['returncode']==0,name
                inv=json.loads((root/name/'coherence.inventory.json').read_text())
                assert inv['sites']==family['retained_sites'],name
                current={tuple(g) for g in inv['groups']}
                if order==1:
                    groups=current
                    assert len(groups)*size==inv['sites']
                    assert sorted(i for g in groups for i in g)==list(range(inv['sites']))
                else:assert current<=groups,('partition changed',name)
                inv.update(case=name,family=family['name'],order=order,sK=size,execution=execution,
                    partition_sha256=hashlib.sha256(json.dumps(sorted(groups)).encode()).hexdigest(),
                    config=cfg,groups_recovered=len(current))
                result['cases'].append(inv);save(root/'summary.json',result)
                print('INVENTORY',name,'clusters',inv['clusters'],'D',inv['max_dimension'],'wall',round(execution['wall_s'],3),flush=True)
    result['status']='passed';save(root/'summary.json',result)


if __name__=='__main__':main()
