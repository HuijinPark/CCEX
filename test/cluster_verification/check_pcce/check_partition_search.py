#!/usr/bin/env python3
"""Build differential assignment checks against the retained full-distance routines."""
import argparse
from pathlib import Path
import subprocess
import unit_build


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
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
        build([cxx,*flags,'test/cluster_verification/check_pcce/test_partition_search.cpp',*objects,*unit_build.link_flags(),'-o',str(root/'checks')])
    with (root/'run.log').open('w') as log:
        subprocess.run([str(root/'checks')],stdout=log,stderr=subprocess.STDOUT,check=True)
    print((root/'run.log').read_text(),end='')


if __name__=='__main__':main()
