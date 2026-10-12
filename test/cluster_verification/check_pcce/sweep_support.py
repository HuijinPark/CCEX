"""Measured, resumable executions for the pCCE factorial experiment."""
import collections
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import threading
import time
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'gpu_verification/check_gpu'))
from benchmark import REPO,compare,parse_profiles
from check_baths import save,sha
from inventory_sweep import resolve


def cluster_data(log,size):
    lines=sorted(line.strip() for line in log.splitlines() if re.match(r'\s*Cluster\[\d+\]\[',line))
    groups=set();by_order=collections.Counter()
    for line in lines:
        m=re.search(r'Cluster\[(\d+)\]\[\s*\d+\]\s*:\s*(-?\d+)\s*\[([^]]+)\]',line)
        if not m:continue
        n,power=int(m[1]),int(m[2]);indices=list(map(int,re.findall(r'\d+',m[3])))
        assert n==len(indices) and n%size==0
        by_order[n//size]+=1
        for j in range(0,n,size):groups.add(tuple(sorted(indices[j:j+size])))
    return dict(cluster_listing_sha256=hashlib.sha256('\n'.join(lines).encode()).hexdigest(),
                groups=sorted(groups),counts_by_order=dict(by_order))


def cgroup_memory():
    """Locate this job step's memory cgroup; report sampled RSS and total usage."""
    for line in Path('/proc/self/cgroup').read_text().splitlines():
        _,controllers,relative=line.split(':',2)
        if 'memory' in controllers.split(','):
            root=Path('/sys/fs/cgroup/memory')/relative.lstrip('/')
            return root,1
        if controllers=='':
            root=Path('/sys/fs/cgroup')/relative.lstrip('/')
            if (root/'memory.current').exists():return root,2
    return None,None


def execute(config,directory,binary,ranks,timeout=1800,monitor=True):
    directory=Path(directory).resolve();binary=Path(binary).resolve()
    config=dict(config,outfile=str(directory/'coherence'))
    existing=directory/'execution.json'
    if existing.exists():
        result=json.loads(existing.read_text())
        assert json.loads((directory/'input.json').read_text())==config,'Resume input changed'
        assert result['binary_sha256']==sha(binary),'Resume binary changed'
        return result
    if directory.exists():
        # Preserve an interrupted Slurm attempt rather than overwriting partial output.
        directory.rename(directory.with_name(directory.name+'.interrupted-'+str(time.time_ns())))
    directory.mkdir(parents=True)
    save(directory/'input.json',config)
    binary_hash=sha(binary)
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    env.setdefault('I_MPI_FABRICS','shm')
    command=['/usr/bin/time','-v','-o',str(directory/'resource.txt'),'mpirun','-n',str(ranks),str(binary),'-f',str(directory/'input.json'),'-v']
    gpu=None;gpu_out=None
    if monitor:
        gpu_out=(directory/'gpu-memory.csv').open('w')
        gpu=subprocess.Popen(['nvidia-smi','--query-gpu=memory.used','--format=csv,noheader,nounits','--loop-ms=100'],stdout=gpu_out,stderr=subprocess.DEVNULL,start_new_session=True)
        # Start monitoring outside the measured interval, including its first query.
        deadline=time.monotonic()+5
        while (directory/'gpu-memory.csv').stat().st_size==0 and gpu.poll() is None and time.monotonic()<deadline:time.sleep(.02)
    cg,version=cgroup_memory() if monitor else (None,None)
    peaks=dict(cgroup_rss_bytes=0,cgroup_usage_bytes=0,samples=0)
    stop=threading.Event()
    def sample():
        while not stop.is_set():
            if cg:
                try:
                    stats=dict(line.split() for line in (cg/'memory.stat').read_text().splitlines())
                    rss=int(stats.get('total_rss',stats.get('rss','0'))) if version==1 else int(stats.get('anon','0'))
                    usage=int((cg/('memory.usage_in_bytes' if version==1 else 'memory.current')).read_text())
                    peaks['cgroup_rss_bytes']=max(peaks['cgroup_rss_bytes'],rss)
                    peaks['cgroup_usage_bytes']=max(peaks['cgroup_usage_bytes'],usage);peaks['samples']+=1
                except (OSError,ValueError):pass
            stop.wait(.1)
    sampler=threading.Thread(target=sample,daemon=True);sampler.start()
    expired=threading.Event()
    with (directory/'run.log').open('w') as log:
        start=time.perf_counter()
        process=subprocess.Popen(command,cwd=REPO,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        def terminate():
            if process.poll() is None:
                expired.set()
                try:os.killpg(process.pid,signal.SIGTERM)
                except ProcessLookupError:return
                # Open MPI normally terminates its ranks on SIGTERM.
                time.sleep(3)
                try:os.killpg(process.pid,signal.SIGKILL)
                except ProcessLookupError:pass
        watchdog=threading.Timer(timeout,terminate);watchdog.daemon=True;watchdog.start()
        code=process.wait();wall=time.perf_counter()-start;watchdog.cancel()
    stop.set();sampler.join()
    if gpu:
        gpu.terminate()
        try:gpu.wait(timeout=5)
        except subprocess.TimeoutExpired:gpu.kill();gpu.wait()
        gpu_out.close()
    log=(directory/'run.log').read_text();profiles=parse_profiles(log)
    counts=collections.Counter()
    for item in profiles:counts[item['dim']]+=int(item['clusters'])
    resource=(directory/'resource.txt').read_text() if (directory/'resource.txt').exists() else ''
    rss=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',resource)
    samples=[]
    if monitor:
        samples=[int(line.strip()) for line in (directory/'gpu-memory.csv').read_text().splitlines() if line.strip().isdigit()]
    clustering=None
    if '       Clusterize ' in log and '       Main calculation ' in log:
        section=log.split('       Clusterize ',1)[1].split('       Main calculation ',1)[0]
        times=re.findall(r'Wall time = ([0-9.eE+-]+) s',section)
        if times:clustering=float(times[-1])
    result=dict(status='timeout' if expired.is_set() else 'passed' if code==0 and profiles else 'failed',
        directory=str(directory),wall_s=wall,returncode=code,binary_sha256=binary_hash,command=command,
        profiles=profiles,cluster_counts_by_dimension=dict(counts),clustering_including_partition_s=clustering,
        memory=dict(**peaks,cgroup=str(cg) if cg else None,cgroup_version=version,
            gnu_time_maxrss_kib=int(rss[1]) if rss else None,gpu_sampled_peak_mib=max(samples) if samples else None,
            gpu_baseline_mib=samples[0] if samples else None,gpu_samples=len(samples),sampling_interval_s=.1,
            note='cgroup RSS/anon and total usage include the Python runner/job step; total usage includes file cache. GNU time maxrss is not the summed simultaneous MPI RSS. GPU samples are whole-device usage.'))
    result.update(cluster_data(log,config['sK']))
    result['raw_sha256']={p.name:sha(p) for suffix in ['*_noDiv','*_wiDiv'] for p in directory.glob(suffix)}
    save(existing,result)
    return result


def check(execution,inventory):
    assert execution['status']=='passed',execution['directory']
    expected=collections.Counter()
    for row in inventory['counts']:expected[str(row['dimension'])]+=row['clusters']
    assert execution['cluster_counts_by_dimension']==dict(expected),'Inventory/profile count mismatch'
    assert execution['cluster_listing_sha256']==inventory['execution']['cluster_listing_sha256'],'Cluster listing changed'
    assert set(map(tuple,execution['groups']))==set(map(tuple,inventory['groups'])),'Partition changed'


def config_for(family,fixture,order,size,pulses,steps=129,seed=None,state=0):
    return dict(resolve(family['config'],fixture),order=order,sK=size,npulse=pulses,nstep=steps,
                deltat=family['end_ms']/(steps-1),seed=family['config']['seed'] if seed is None else seed,nstate=state)


def measure(config,directory,ranks,repeats,binaries,inventory=None):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    result=dict(status='running',config=config,ranks=ranks,trials=[])
    for repeat in range(repeats):
        variants=['cpu','cuda','auto'];variants=variants[repeat%3:]+variants[:repeat%3]
        trial=dict(repeat=repeat,executions={},comparisons={});result['trials'].append(trial)
        for variant in variants:
            execution=execute(dict(config,backend=variant),directory/f'repeat{repeat}'/variant,
                              binaries['cpu' if variant=='cpu' else 'cuda'],ranks)
            trial['executions'][variant]=execution;save(directory/'summary.json',result)
            if execution['status']!='passed':
                result.update(status='incomplete',reason=execution['status']);save(directory/'summary.json',result);return result
            if inventory:check(execution,inventory)
        cpu=trial['executions']['cpu']
        for variant,execution in trial['executions'].items():
            for key in ['cluster_listing_sha256','cluster_counts_by_dimension','groups']:
                assert cpu[key]==execution[key],key
            trial['comparisons'][variant]=compare(Path(cpu['directory']),Path(execution['directory']),1e-8,2e-6)
        print('PASS',directory.name,'MPI',ranks,'repeat',repeat,{v:round(e['wall_s'],3) for v,e in trial['executions'].items()},flush=True)
    import statistics
    result['median_s']={v:statistics.median(t['executions'][v]['wall_s'] for t in result['trials']) for v in ['cpu','cuda','auto']}
    result['status']='passed';save(directory/'summary.json',result)
    return result
