"""How the unit checks here are built: against the objects of the current build, with its compiler and flags,
as the build itself gives them through test/unit_build_env.sh (D26, 2026-10-02): the Makefile of do_compile.sh
(its unit-env target), or a CMake build of any generator (CCEX_BUILD_DIR=<its build folder>, its unit-build.env).
That script brings the build up to date first, so the objects are those of the current sources.
The tests add -UNDEBUG for their asserts, and --sanitize its own flags. Until 2026-10-02 this module read the
Makefile, CMakeCache.txt and flags.make itself (Unix Makefiles generator only) and compared every object with its
sources and headers; before that it added -march=native whatever the build had (D19)."""
import os
import shlex
import subprocess
from functools import lru_cache
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
SOURCES = {'cluster_pcce': 'src/clustering/cluster_pcce.cpp', 'cluster_hash': 'src/clustering/cluster_hash.cpp',
           'cluster': 'src/clustering/cluster.cpp', 'memory': 'src/common/memory.cpp'}


@lru_cache(maxsize=None)
def _env():
    out = subprocess.run(['bash', str(REPO/'test'/'unit_build_env.sh')], capture_output=True, text=True)
    if out.returncode:
        raise SystemExit(out.stderr.strip() or 'test/unit_build_env.sh failed')
    values = {}
    for line in out.stdout.splitlines():
        key, _, value = line.partition('=')
        if key.startswith('CCEX_UNIT_'):
            values[key] = value.strip("'")
    return values


def compiler():
    # what the compiler needs in its environment (I_MPI_CXX=icpx on idm): the checks call it as a subprocess
    for pair in _env()['CCEX_UNIT_ENV'].split():
        name, _, value = pair.partition('=')
        os.environ[name] = value
    return _env()['CCEX_UNIT_CXX']


def objects():
    found = _env()['CCEX_UNIT_OBJS'].split()
    if not found:
        raise SystemExit('test/unit_build_env.sh gave no objects')
    return found


def without(found, name):
    """found less the object of source `name` (e.g. 'cluster_pcce'), whichever build it comes from"""
    return [p for p in found if Path(p).name not in (f'{name}.o', f'{name}.cpp.o')]


def build_flags():
    """the compile flags and includes of the build whose objects are linked"""
    return shlex.split(_env()['CCEX_UNIT_CXXFLAGS']) + shlex.split(_env()['CCEX_UNIT_INCLUDES'])


def flags(sanitize):
    base = build_flags()+['-UNDEBUG', '-fpermissive', '-Wno-deprecated-declarations',
                          '-Iinclude', '-Iexternal/eigen', '-Iexternal/uthash/include']
    return base + (['-O1', '-fsanitize=address,undefined', '-fno-omit-frame-pointer'] if sanitize else [])


def link_flags():
    return shlex.split(_env()['CCEX_UNIT_LIBS'])
