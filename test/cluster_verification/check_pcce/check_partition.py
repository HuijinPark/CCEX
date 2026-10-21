#!/usr/bin/env python3
"""Build/run partition unit checks against the current do_compile.sh build.

--sanitize instruments the partition, hash, cluster ownership and allocation
code with AddressSanitizer/UndefinedBehaviorSanitizer; leak detection is enabled
by the sanitizer runtime. Results are written to a new output directory.
"""
import argparse
from pathlib import Path
import subprocess
import unit_build


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mpicxx', default=None, help='the compiler; default: the Makefile\'s CXX')
    parser.add_argument('--sanitize', action='store_true')
    args = parser.parse_args()
    repo = unit_build.REPO
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    objects = unit_build.objects()
    flags = unit_build.flags(args.sanitize)
    cxx = args.mpicxx or unit_build.compiler()
    with (output/'build.log').open('w') as log:
        def build(command):
            subprocess.run(command, cwd=repo, stdout=log, stderr=subprocess.STDOUT, check=True)
        if args.sanitize:
            flags += ['-ffunction-sections', '-fdata-sections']
            for name, source in unit_build.SOURCES.items():
                obj = output/f'{name}.o'
                build([cxx, *flags, '-c', source, '-o', str(obj)])
                objects = unit_build.without(objects, name) + [str(obj)]
        binary = output/'partition-checks'
        build([cxx, *flags, 'test/cluster_verification/check_pcce/test_partition_checks.cpp', *objects,
               *unit_build.link_flags(), '-Wl,--gc-sections', '-o', str(binary)])
    with (output/'run.log').open('w') as log:
        subprocess.run([str(binary)], stdout=log, stderr=subprocess.STDOUT, check=True)
    with (output/'initialization.txt').open('w') as log:
        subprocess.run([str(binary), 'snapshot'], stdout=log, check=True)
    print(f'PASS partition checks: {output}')


if __name__ == '__main__':
    main()
