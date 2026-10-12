#!/usr/bin/env python3
"""Exercise the production pCCE cache, including MPI, sampled states and corrupt files."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import struct
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'gpu_verification/check_gpu'))
from check_pcce import run_case
from benchmark import REPO,compare
from partition_log import check_search,check_cache_hit


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repair_checksum(data):
    result=14695981039346656037
    for byte in data[:-8]:
        result=((result ^ byte)*1099511628211)&((1<<64)-1)
    data[-8:]=struct.pack('<Q',result)
    return data


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cpu',type=Path,required=True)
    p.add_argument('--cuda',type=Path)
    p.add_argument('--baseline',type=Path)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    binaries={'cpu':a.cpu.resolve()}
    if a.cuda:binaries['cuda']=a.cuda.resolve()
    if a.baseline:binaries['baseline']=a.baseline.resolve()
    original=json.loads((REPO/'test/gpu_verification/check_gpu/examples/pcce_ensemble.json').read_text())
    original['gyrofile']=str(REPO/original['gyrofile'])
    original['bathfile']=[str(REPO/name) for name in original['bathfile']]
    original.update(backend='cpu',iter_detail=True,max_trial=5,pcce_partition='deterministic')
    records=[]
    summary=dict(status='running',binary_sha256={k:sha(v) for k,v in binaries.items()},records=records)

    def save():
        (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')

    def run(name,config,binary='cpu',ranks=1,event=None,reference=None,error=None,extra=()):
        directory=root/name
        execution,log=run_case(config,directory,binaries[binary],ranks,180,extra=extra)
        record=dict(case=name,ranks=ranks,binary=binary,execution=execution)
        if error:
            assert execution['returncode']!=0 and error in log,(name,log[-3000:])
            record['expected_error']=error
        else:
            assert execution['returncode']==0,(name,log[-3000:])
            if event:
                assert f'[partition-cache] status={event}' in log,(name,event)
            if event=='hit':
                check_cache_hit(log)
            elif config.get('pcce_partition')!='legacy':
                record['trial_audit']=check_search(log,config['max_trial'],ranks,config.get('pcce_detail_limit',32))
            if reference:
                ref=json.loads((reference/'execution.json').read_text())
                for key in ['cluster_listing_sha256','cluster_counts_by_dimension']:
                    assert execution[key]==ref[key],(name,key)
                record['comparison']=compare(reference,directory,1e-8,2e-6)
            else:record['comparison']=compare(directory,directory,1e-8,2e-6)
        records.append(record);save();print('PASS',name,flush=True)
        return directory

    try:
        # Reuse across order, rdip, pulses, ensemble/sampled states and ranks;
        # these change subsequent dynamics/cluster construction, not partition.
        for size in [1,2,3]:
            base=dict(original,sK=size)
            cache=root/f'partition-k{size}.cache'
            cold=dict(base,pcce_cache_mode='readwrite',pcce_cache_file=str(cache))
            reference=run(f'k{size}/reference',base)
            if a.baseline:run(f'k{size}/before',base,'baseline',reference=reference)
            run(f'k{size}/create',cold,ranks=2,event='stored',reference=reference)
            checksum=sha(cache)
            variants=[dict(nstate=0,npulse=1,order=2,rdip=30),
                      dict(nstate=2,npulse=4,order=2,rdip=30),
                      dict(nstate=0,npulse=16,order=1,rdip=30),
                      dict(nstate=1,npulse=1,order=2,rdip=7)]
            for index,changes in enumerate(variants):
                cfg=dict(base,**changes)
                ref=reference if index==0 else run(f'k{size}/v{index}/reference',cfg)
                cached=dict(cfg,pcce_cache_mode='read',pcce_cache_file=str(cache))
                for ranks in [1,2,4]:
                    run(f'k{size}/v{index}/cpu-mpi{ranks}',cached,ranks=ranks,event='hit',reference=ref)
                if a.cuda:
                    for backend in ['cuda','auto']:
                        run(f'k{size}/v{index}/{backend}',dict(cached,backend=backend),'cuda',2,'hit',ref)
                assert sha(cache)==checksum,'read mode modified the cache'

        base=dict(original,sK=2)
        clean=(root/'partition-k2.cache').read_bytes()
        mutations=[('seed',{'seed':73}),('sK',{'sK':3}),('trials',{'max_trial':7}),
                   ('iterations',{'max_iter':11}),('initialization',{'kmeans_pp':False})]
        bath=Path(original['bathfile'][0]).read_text().splitlines()
        for label,rows in [('reordered',[bath[0],*reversed(bath[1:])]),
                           ('coordinate',[bath[0],bath[1].replace('5.0','5.125'),*bath[2:]]),
                           ('remainder',['8.0 0 0 0',*bath[1:],'1.0 1.0 1.0 13C'])]:
            path=root/f'bath-{label}';path.write_text('\n'.join(rows)+'\n')
            mutations.append((label,{'bathfile':[str(path)]}))
        for name,changes in mutations:
            cache=root/f'mismatch-{name}.cache';cache.write_bytes(clean)
            cfg=dict(base,**changes)
            ref=run(f'invalidate/{name}/reference',cfg)
            cached=dict(cfg,pcce_cache_file=str(cache),pcce_cache_mode='read')
            run(f'invalidate/{name}/reject',cached,ranks=2,error='read mode requires a valid matching cache')
            assert cache.read_bytes()==clean
            run(f'invalidate/{name}/replace',dict(cached,pcce_cache_mode='readwrite'),ranks=2,event='stored',reference=ref)
            assert cache.read_bytes()!=clean,name
            run(f'invalidate/{name}/hit',cached,ranks=4,event='hit',reference=ref)

        # Format is magic + length-prefixed complete identity + binary64 centers
        # + uint64 assignments + FNV-1a64 corruption checksum, all little endian.
        magic=b'CCEX-PCCE-CACHE-v1\n'
        start=len(magic)+8+struct.unpack_from('<Q',clean,len(magic))[0]
        bad_center=bytearray(clean);struct.pack_into('<Q',bad_center,start,0x7ff8000000000000)
        bad_index=bytearray(clean);struct.pack_into('<Q',bad_index,start+3*3*8,3)
        bad_size=bytearray(clean)
        assignments=start+3*3*8
        first=struct.unpack_from('<Q',clean,assignments)[0]
        struct.pack_into('<Q',bad_size,assignments,(first+1)%3)
        flipped=bytearray(clean);flipped[-10]^=1
        wrong_version=bytearray(clean);wrong_version[len(magic)+16]^=1
        corruptions={'truncated':clean[:40],'oversized':clean+b'junk','bitflip':flipped,
                     'nonfinite-center':repair_checksum(bad_center),
                     'invalid-index':repair_checksum(bad_index),'group-size':repair_checksum(bad_size),
                     'algorithm-version':repair_checksum(wrong_version)}
        ref=root/'k2/reference'
        for name,data in corruptions.items():
            cache=root/f'corrupt-{name}.cache';cache.write_bytes(data)
            cfg=dict(base,pcce_cache_file=str(cache),pcce_cache_mode='read')
            run(f'corrupt/{name}/reject',cfg,ranks=2,error='read mode requires a valid matching cache')
            assert cache.read_bytes()==data
            run(f'corrupt/{name}/repair',dict(cfg,pcce_cache_mode='readwrite'),ranks=2,event='stored',reference=ref)
            assert cache.read_bytes()==clean,name

        cache=root/'partition-k2.cache'
        active=dict(base,pcce_cache_mode='read',pcce_cache_file=str(cache))
        for name,changes,message,extra in [
            ('missing',{'pcce_cache_file':str(root/'missing.cache')},'read mode requires a valid matching cache',()),
            ('path-required',{'pcce_cache_file':None},'enabled cache requires pcce_cache_file',()),
            ('mode-type',{'pcce_cache_mode':True},'pcce_cache_mode must be',()),
            ('bad-mode',{'pcce_cache_mode':'write'},'pcce_cache_mode must be',()),
            ('path-type',{'pcce_cache_file':3},'pcce_cache_file must be',()),
            ('empty-path',{'pcce_cache_file':''},'pcce_cache_file must be',()),
            ('long-path',{'pcce_cache_file':'a'*500},'pcce_cache_file must be',()),
            ('legacy',{'pcce_partition':'legacy'},'enabled cache requires pcce_partition=deterministic',()),
            ('method',{'method':'CCE'},'enabled cache requires method=pCCE',()),
            ('cli-method',{},'enabled cache requires method=pCCE',('-m','CCE')),
            ('directory',{'pcce_cache_file':str(root)},'cache path must be a regular file',()),
            ('write-failure',{'pcce_cache_mode':'readwrite','pcce_cache_file':str(root/'absent-dir/p.cache')},'cannot create cache temporary file',()),
        ]:
            cfg=dict(active,**changes)
            cfg={k:v for k,v in cfg.items() if v is not None}
            run(f'options/{name}',cfg,ranks=2,error=message,extra=extra)

        # Explicit off still has the historical RNG behavior, including legacy.
        checksum=sha(cache)
        run('off-ignores-path',dict(active,pcce_cache_mode='off'),reference=ref)
        assert sha(cache)==checksum
        if a.baseline:
            # Compare the historical seed-zero C RNG stream, not the old bug
            # that ignored every explicitly requested seed.
            cfg=dict(base,pcce_partition='legacy',nstate=2,seed=0)
            ref_legacy=run('legacy/before',cfg,'baseline',ranks=2)
            run('legacy/after',dict(cfg,pcce_cache_mode='off',pcce_cache_file=str(cache)),ranks=2,reference=ref_legacy)

        # Two independent jobs can publish the same file without exposing a
        # partially written partition. Compare every raw result with uncached.
        concurrent=root/'concurrent.cache'
        cfg=dict(base,pcce_cache_mode='readwrite',pcce_cache_file=str(concurrent))
        def writer(index):
            directory=root/f'concurrent/writer{index}'
            execution,log=run_case(cfg,directory,binaries['cpu'],1,180)
            assert execution['returncode']==0,log[-2000:]
            return dict(case=f'concurrent/writer{index}',execution=execution,comparison=compare(ref,directory,1e-8,2e-6))
        with ThreadPoolExecutor(max_workers=2) as pool:
            records.extend(pool.map(writer,range(2)))
        assert concurrent.read_bytes()==clean
        run('concurrent/read',dict(cfg,pcce_cache_mode='read'),ranks=4,event='hit',reference=ref)
        assert not list(root.glob('*.tmp.*')),'temporary cache file leaked'
        summary['status']='passed';save()
    except BaseException as error:
        summary.update(status='failed',error=f'{type(error).__name__}: {error}');save();raise


if __name__=='__main__':main()
