#!/usr/bin/env python3
"""Check cycling real baths, MPI statistics, legacy RNG, and cache migration."""
import argparse
import json
import math
from pathlib import Path
import re
import struct
import sys
from sweep_support import REPO,config_for,compare,sha,save
sys.path.insert(0,str(REPO/'test/gpu_verification/check_gpu'))
from check_pcce import run_case


def stats(log):
    match=re.search(r'\[partition-search\] converged=(\d+) cycles=(\d+) limit=(\d+) iterations=(\d+) executed=(\d+) skipped=(\d+)',log)
    assert match,'Missing search statistics'
    result=dict(zip(['converged','cycles','limit','iterations','executed','skipped'],map(int,match.groups())))
    assert result['iterations']==result['executed']+result['skipped']
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fixtures',type=Path,required=True)
    p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    fixture=a.fixtures.resolve();family=json.loads((fixture/'manifest.json').read_text())['families'][0]
    cfg=config_for(family,fixture,2,3,4,steps=9,state=1)
    cfg.update(max_trial=7,max_iter=1001,iter_detail=True,pcce_cache_mode='off')
    summary=dict(status='running',records=[])
    def run(name,config,ranks=2,backend='cpu',baseline=False,reference=None,error=None,cache_hit=False):
        binary=a.baseline.resolve() if baseline else REPO/f"bin/main.{'cpu' if backend=='cpu' else 'cuda'}.out"
        directory=root/name
        execution,log=run_case(dict(config,backend=backend),directory,binary,ranks,180)
        record=dict(case=name,ranks=ranks,backend=backend,baseline=baseline,execution=execution,binary_sha256=sha(binary))
        if error:
            assert execution['returncode']!=0 and error in log,(name,log[-3000:])
            record['expected_error']=error
        else:
            assert execution['returncode']==0,(name,log[-3000:])
            if not baseline and not cache_hit:
                record['search']=stats(log)
                count=config['max_trial'] if config.get('pcce_partition')!='legacy' else max(1,config['max_trial']//ranks)*ranks
                assert sum(record['search'][key] for key in ['converged','cycles','limit'])==count
            if reference:
                ref=json.loads((reference/'execution.json').read_text())
                for key in ['cluster_listing_sha256','cluster_counts_by_dimension']:assert ref[key]==execution[key],(name,key)
                record['comparison']=compare(reference,directory,1e-8,2e-6)
        summary['records'].append(record);save(root/'summary.json',summary)
        print('PASS',name,record.get('search',{}),flush=True)
        return directory,log
    try:
        reference,_=run('nv-before',cfg,ranks=1,baseline=True)
        reference_stats=None
        for ranks,backend in [(1,'cpu'),(2,'cpu'),(4,'cpu'),(2,'cuda'),(2,'auto')]:
            _,log=run(f'nv-{backend}-mpi{ranks}',cfg,ranks,backend,reference=reference)
            current=stats(log);assert current['cycles']>0 and current['skipped']>0
            if reference_stats is None:reference_stats=current
            else:assert current==reference_stats,'MPI/backend changed deterministic search work'
        legacy=dict(cfg,pcce_partition='legacy',max_trial=8,max_iter=101,nstate=2)
        ref,_=run('legacy-before',legacy,baseline=True)
        for backend in ['cpu','cuda']:
            _,log=run(f'legacy-{backend}',legacy,backend=backend,reference=ref)
            assert stats(log)['cycles']>0 and stats(log)['skipped']>0

        idle=dict(cfg,max_trial=1)
        ref,_=run('idle-before',idle,ranks=1,baseline=True)
        run('idle-mpi4',idle,ranks=4,reference=ref)

        # Observe short periods first, then independently run the old loop to the
        # phase that an INT_MAX run reaches. Never run billions of old iterations.
        probe=dict(cfg,max_trial=3,max_iter=101)
        _,log=run('large-counter-probe',probe,ranks=2)
        periods=list(map(int,re.findall(r'\[partition-iteration\][^\n]*status=cycle[^\n]*period=(\d+)',log)))
        assert len(periods)>=2
        common=math.lcm(*periods);assert common<10000
        end=2147483647
        phase=101+(end-101)%common
        ref,_=run('large-counter-before-phase',dict(probe,max_iter=phase),ranks=2,baseline=True)
        _,log=run('large-counter-mpi4',dict(probe,max_iter=end),ranks=4,reference=ref)
        assert stats(log)['iterations']>4294967295 and stats(log)['executed']<1000

        # Old and new files carry different algorithm identities but must have
        # identical exact centers/assignments after their freshly computed search.
        cache=root/'partition.cache'
        cold=dict(cfg,pcce_cache_mode='readwrite',pcce_cache_file=str(cache),nstate=0)
        ref,_=run('cache-v2-before',cold,baseline=True)
        old=cache.read_bytes();assert b'pcce-deterministic-20260912-v2' in old
        run('cache-v2-read-rejected',dict(cold,pcce_cache_mode='read'),error='read mode requires a valid matching cache')
        assert cache.read_bytes()==old
        run('cache-v3-replace',cold,reference=ref)
        new=cache.read_bytes();assert b'pcce-deterministic-20260912-v3' in new
        magic=b'CCEX-PCCE-CACHE-v1\n'
        def payload(data):return data[len(magic)+8+struct.unpack_from('<Q',data,len(magic))[0]:-8]
        assert payload(old)==payload(new)
        (root/'partition-v2.cache').write_bytes(old)
        run('cache-v3-gpu-read',dict(cold,pcce_cache_mode='read'),backend='cuda',reference=ref,cache_hit=True)
        assert cache.read_bytes()==new
        summary.update(status='passed',exact_cache_payload_equal=True)
        save(root/'summary.json',summary)
    except BaseException as error:
        summary.update(status='failed',error=str(error));save(root/'summary.json',summary);raise


if __name__=='__main__':main()
