#!/usr/bin/env python3
"""Compare cold partition search and verified reuse, rotating CPU/CUDA/auto order."""
import argparse
import json
from pathlib import Path
import statistics
from sweep_support import REPO,execute,config_for,compare,sha,save
from benchmark import put


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fixtures',type=Path,required=True)
    p.add_argument('--inventory',type=Path,required=True)
    p.add_argument('--reference-root',type=Path,help='Unmodified pilot result directory for raw comparisons')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--repeats',type=int,default=3)
    a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    fixture=a.fixtures.resolve();manifest=json.loads((fixture/'manifest.json').read_text())
    for name,want in manifest['files'].items():assert sha(fixture/name)==want,name
    inventory=json.loads(a.inventory.read_text());assert inventory['status']=='passed'
    inventories={case['case']:case for case in inventory['cases']}
    summary=dict(status='running',ranks=4,repeats=a.repeats,cases=[],
                 manifest_sha256=sha(fixture/'manifest.json'),binary_sha256={backend:sha(REPO/f'bin/main.{backend}.out') for backend in ['cpu','cuda']})
    save(root/'summary.json',summary)
    try:
        for family in manifest['families']:
            # NV includes D512 GPU work; PL6 isolates the measured costly partition.
            order,size=(3,3) if family['name']=='nv13c' else (2,2)
            name=f"{family['name']}-o{order}-k{size}-p1"
            inv=inventories[f"{family['name']}-o{order}-k{size}"]
            config=config_for(family,fixture,order,size,1,steps=33)
            case=dict(case=name,config=config,sites=family['retained_sites'],species=family['species'],
                      inventory=inv,executions=[],medians={})
            summary['cases'].append(case)
            reference=None;cached_file=None
            from sweep_support import check
            for phase in ['cold','reuse']:
                for repeat in range(a.repeats):
                    backends=['cpu','cuda','auto'];backends=backends[repeat%3:]+backends[:repeat%3]
                    for backend in backends:
                        directory=root/name/f'{phase}-{repeat}-{backend}'
                        cache=root/name/f'partition-{repeat}-{backend}.cache' if phase=='cold' else cached_file
                        cfg=put(config,{'backend':backend,'pcce.cache.mode':'readwrite' if phase=='cold' else 'read','pcce.cache.file':str(cache)})
                        binary=REPO/f"bin/main.{'cpu' if backend=='cpu' else 'cuda'}.out"
                        before=sha(cache) if phase=='reuse' else None
                        run=execute(cfg,directory,binary,4,timeout=1800)
                        assert run['status']=='passed',(name,phase,backend,run)
                        check(run,inv)
                        log=(directory/'run.log').read_text()
                        assert f"[partition-cache] status={'stored' if phase=='cold' else 'hit'}" in log
                        if reference is None:
                            reference=directory;cached_file=cache
                            if a.reference_root:
                                old=a.reference_root/name/'repeat0/cpu'
                                assert old.exists(),old
                                case['prior_pilot_comparison']=compare(old,directory,1e-8,2e-6)
                        error=compare(reference,directory,1e-8,2e-6)
                        if phase=='reuse':assert sha(cache)==before
                        case['executions'].append(dict(phase=phase,repeat=repeat,backend=backend,execution=run,comparison=error))
                        save(root/'summary.json',summary)
                        print('PASS',name,phase,backend,repeat,'wall',round(run['wall_s'],4),flush=True)
            for backend in ['cpu','cuda','auto']:
                times={phase:statistics.median(r['execution']['wall_s'] for r in case['executions'] if r['phase']==phase and r['backend']==backend) for phase in ['cold','reuse']}
                times['cold_over_reuse']=times['cold']/times['reuse'];case['medians'][backend]=times
            save(root/'summary.json',summary)
        summary['status']='passed';save(root/'summary.json',summary)
    except BaseException as error:
        summary.update(status='failed',error=f'{type(error).__name__}: {error}');save(root/'summary.json',summary);raise


if __name__=='__main__':main()
