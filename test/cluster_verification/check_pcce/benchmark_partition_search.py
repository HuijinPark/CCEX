#!/usr/bin/env python3
"""Compare fresh pCCE searches on different physical bath realizations; cache off."""
import argparse
import json
from pathlib import Path
import re
import statistics
import struct
from sweep_support import REPO,execute,config_for,compare,sha,save
from benchmark import get,put


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fixtures',type=Path,required=True)
    p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ranks',type=int,default=32)
    p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--scope',choices=['full','realizations-cpu'],default='full')
    a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    fixture=a.fixtures.resolve();baseline=a.baseline.resolve()
    manifest=json.loads((fixture/'search-manifest.json').read_text())
    for name,want in manifest['files'].items():assert sha(fixture/name)==want,name
    for item in manifest['realizations']:assert sha(fixture/item['name'])==item['sha256']
    families={f['name']:f for f in manifest['families']}
    summary=dict(status='running',scope=a.scope,note='Every timed execution has pcce.cache.mode off; no partition reuse.',
                 manifest=manifest,cases=[],exact_partitions=[],binary_sha256={})
    for version,base in [('before',baseline),('after',REPO/'bin')]:
        for backend in ['cpu','cuda']:
            summary['binary_sha256'][version+'-'+backend]=sha(base/f'main.{backend}.out')
    save(root/'summary.json',summary)

    def run_case(name,cfg,ranks,repeats,backends):
        case=dict(case=name,config=cfg,ranks=ranks,repeats=repeats,executions=[],medians={})
        summary['cases'].append(case);save(root/'summary.json',summary)
        reference=None
        for repeat in range(repeats):
            for backend in backends:
                for version in (['before','after'] if repeat%2==0 else ['after','before']):
                    binary=(baseline if version=='before' else REPO/'bin')/f"main.{'cpu' if backend=='cpu' else 'cuda'}.out"
                    directory=root/name/f'{backend}-{repeat}-{version}'
                    run=execute(put(cfg,{'backend':backend,'pcce.cache.mode':'off'}),directory,binary,ranks,timeout=1200)
                    assert run['status']=='passed',(name,backend,version,run)
                    log=(directory/'run.log').read_text()
                    assert '[partition-cache]' not in log and 'source=cache' not in log
                    trial=re.search(r'\[partition\] mode=deterministic total_trials=(\d+) best_trial=(\d+)',log)
                    assert trial and int(trial[1])==get(cfg,'pcce.trials')
                    run['best_trial']=int(trial[2])
                    if reference is None:reference=(directory,run)
                    for key in ['groups','cluster_listing_sha256','cluster_counts_by_dimension','best_trial']:
                        assert reference[1][key]==run[key],(name,version,key)
                    error=compare(reference[0],directory,1e-8,2e-6)
                    case['executions'].append(dict(backend=backend,version=version,repeat=repeat,execution=run,comparison=error))
                    save(root/'summary.json',summary)
                    print('PASS',name,backend,version,repeat,'wall',round(run['wall_s'],4),flush=True)
        for backend in backends:
            times={v:statistics.median(e['execution']['wall_s'] for e in case['executions'] if e['backend']==backend and e['version']==v) for v in ['before','after']}
            times['speedup']=times['before']/times['after'];case['medians'][backend]=times
        save(root/'summary.json',summary)

    try:
        # Actual supplied bath realizations, with each realization's complete
        # radius-selected bath. Standard pCCE remainder removal remains enabled.
        for item in manifest['realizations']:
            cfg=config_for(families['pl6'],fixture,2,2,1,steps=33)
            cfg['configuration']['bathfile']=[str(fixture/item['name'])]
            backends=['cpu','cuda','auto'] if item['realization']==1 and a.scope=='full' else ['cpu']
            run_case(item['name']+'-k2',cfg,a.ranks,a.repeats,backends)
        if a.scope=='realizations-cpu':
            summary['status']='passed';save(root/'summary.json',summary);return
        # Group sizes and a second material check that the optimization is general.
        for family,sizes in [('nv13c',[1,2,3]),('pl6',[1,3])]:
            for size in sizes:
                cfg=config_for(families[family],fixture,2,size,4,steps=33)
                run_case(f'{family}-k{size}-pulse4',cfg,a.ranks,1,['cpu','cuda'])
        cfg=config_for(families['pl6'],fixture,2,2,1,steps=33)
        run_case('pl6-fixed786-mpi4',cfg,4,1,['cpu'])

        # Verify binary64 center coordinates as well as group IDs using freshly
        # generated caches. Version identity intentionally differs after the edit.
        for family,size in [('nv13c',3),('pl6',2)]:
            cfg=config_for(families[family],fixture,2,size,1,steps=3)
            payloads=[]
            for version in ['before','after']:
                cache=root/f'{family}-{version}.partition'
                directory=root/f'exact-{family}-{version}'
                binary=(baseline if version=='before' else REPO/'bin')/'main.cpu.out'
                run=execute(put(cfg,{'pcce.cache.mode':'readwrite','pcce.cache.file':str(cache)}),directory,binary,a.ranks)
                assert run['status']=='passed'
                data=cache.read_bytes();magic=b'CCEX-PCCE-CACHE-v1\n'
                assert data.startswith(magic)
                offset=len(magic)+8+struct.unpack_from('<Q',data,len(magic))[0]
                payloads.append(data[offset:-8])
            assert payloads[0]==payloads[1],family
            summary['exact_partitions'].append(dict(family=family,sK=size,binary64_centers_and_assignments_equal=True))
            save(root/'summary.json',summary)
        summary['status']='passed';save(root/'summary.json',summary)
    except BaseException as error:
        summary.update(status='failed',error=f'{type(error).__name__}: {error}')
        save(root/'summary.json',summary);raise


if __name__=='__main__':main()
