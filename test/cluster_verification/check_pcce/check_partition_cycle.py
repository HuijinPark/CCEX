#!/usr/bin/env python3
"""Verify exact final pCCE partitions with and without cycle skipping."""
import argparse
import json
from pathlib import Path
import subprocess
import unit_build


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--fixtures',type=Path)
    p.add_argument('--sanitize',action='store_true')
    a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    repo=unit_build.REPO
    objects=unit_build.objects()
    flags=unit_build.flags(a.sanitize)
    cxx=unit_build.compiler()
    with (root/'build.log').open('w') as log:
        def build(cmd):subprocess.run(cmd,cwd=repo,stdout=log,stderr=subprocess.STDOUT,check=True)
        if a.sanitize:
            obj=root/'cluster_pcce.o'
            build([cxx,*flags,'-c',unit_build.SOURCES['cluster_pcce'],'-o',str(obj)])
            objects=unit_build.without(objects,'cluster_pcce')+[str(obj)]
        binary=root/'checks'
        build([cxx,*flags,'test/cluster_verification/check_pcce/test_partition_cycle.cpp',*objects,*unit_build.link_flags(),'-o',str(binary)])
    cases=[('synthetic',[])]
    if a.fixtures:
        f=a.fixtures.resolve()
        cases += [(f'pl6-bath{i}-k{k}',[str(f/f'pl6-bath-realization-{i}'),str(k),'1788434958','16','1000']) for i,k in [(1,2),(2,2),(3,2),(1,3)]]
        cases += [('nv13c-k3',[str(f/'nv13c-bath'),'3','73','24','40'])]
    summary=dict(status='running',sanitize=a.sanitize,cases=[])
    try:
        for name,args in cases:
            with (root/f'{name}.log').open('w') as log:
                subprocess.run([str(binary),*args],stdout=log,stderr=subprocess.STDOUT,check=True)
            output=(root/f'{name}.log').read_text();print(name,output.strip(),flush=True)
            summary['cases'].append(dict(case=name,output=output))
        summary['status']='passed'
    except BaseException as error:
        summary.update(status='failed',error=str(error));raise
    finally:(root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':main()
