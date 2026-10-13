#!/usr/bin/env python3
"""The refinement policies against their own recorded answers. No Python oracle.

    regress.py WORKDIR                 # run every case, compare with golden/
    regress.py WORKDIR --only NAME     # one case
    regress.py WORKDIR --update        # rewrite golden/ from this run - deliberately

Each case in regression/cases.json is run for real by the C++ driver
(test/cluster_verification/check_adaptive/cpp/bin/ccex-adaptive, which launches bin/main.out for every step), and
what it decided - its exit status and its summary, one per condition where there
are several - has to be what regression/golden/ holds. A decision has to match
exactly; a float may move by 1e-9, which is what was settled on 2026-09-18 for the
comparison against Python that this replaces.

What is measured rather than decided is left out before comparing: elapsed times,
the seconds a policy records beside a D^3 price, and the hash of the binary (every
build has its own). Fault cases need the launcher the recordings used,
test/cluster_verification/check_adaptive/fault_mpiexec.py.

The golden files were first written by this driver at the commit where all 57 of
its runs had just been shown to agree with the Python reference, decision for
decision, so they carry that agreement forward.
"""
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
REPO = TOOLS.parent.parent.parent
SET = TOOLS/'regression'
DRIVER = Path(os.environ.get('CCEX_ADAPTIVE_BIN', TOOLS/'cpp/bin/ccex-adaptive'))
BINARY = Path(os.environ.get('CCEX_BIN', REPO/'bin/main.out'))
MEASURED = {'wall_s', 'search_s', 'seconds_used', 'timing', 'measured_local_seconds',
            'measured_worker_seconds', 'binary_sha256', 'context_sha256'}


def normalized(x, swaps):
    if isinstance(x, dict):              # keys too: input_file_sha256 is keyed by the input's path (from a worktree the
        return {normalized(k, swaps): normalized(v, swaps) for k, v in x.items()      # golden's absolute path differed)
                if k not in MEASURED and not k.endswith('_wall_s')}
    if isinstance(x, list):
        return [normalized(v, swaps) for v in x]
    if isinstance(x, str):
        for old, new in swaps:
            x = x.replace(old, new)
    return x


def timed_out_alike(x):
    """A step the driver ended at its timeout records "timeout": true and a nonzero returncode, which is -15 when the
    SIGTERM ended the launcher and 15 when the launcher caught it and exited with it (Intel mpirun does either; the
    golden answers hold -15). The decision reads the timeout only (policy_merge.cpp), so for the comparison every such
    returncode is the same. Until 2026-10-02 fault-probe-timeout and fault-deadline differed run to run on WSL."""
    if isinstance(x, dict):
        y = {k: timed_out_alike(v) for k, v in x.items()}
        if y.get('timeout') is True and isinstance(y.get('returncode'), int) and y['returncode'] != 0:
            y['returncode'] = 'nonzero after the timeout'
        return y
    if isinstance(x, list):
        return [timed_out_alike(v) for v in x]
    return x


def differences(a, b, at=''):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                yield f'{at}/{k}: only in {"golden" if k in a else "this run"}'
            else:
                yield from differences(a[k], b[k], f'{at}/{k}')
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            yield f'{at}: {len(a)} entries in golden, {len(b)} in this run'
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                yield from differences(x, y, f'{at}[{i}]')
    elif isinstance(a, float) or isinstance(b, float):
        ok = isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool) \
            and not isinstance(b, bool) and (a == b or (math.isnan(a) and math.isnan(b)) or
                                             abs(a-b) <= 1e-9*max(1., abs(a), abs(b)))
        if not ok:
            yield f'{at}: {a!r} in golden, {b!r} in this run'
    elif a != b or type(a) is not type(b):
        yield f'{at}: {a!r} in golden, {b!r} in this run'


def run(case, work):
    out = work/case['name']
    shutil.rmtree(out, ignore_errors=True)
    for stale in work.glob(f"{case['name']}.*"):
        stale.unlink()
    fill = {'<INPUTS>': str(SET/'inputs'), '<WORK>': str(work), '<TOOLS>': str(TOOLS),
            '<BINARY>': str(BINARY), '<REPO>': str(REPO)}

    def filled(text):
        for key, value in fill.items():
            text = text.replace(key, value)
        return text
    env = dict(os.environ, **{k: filled(v) for k, v in case['env'].items()})
    done = subprocess.run([str(DRIVER), *map(filled, case['args'])], env=env, cwd=str(REPO),
                          capture_output=True, text=True)
    swaps = [(str(out), '<OUT>'), (str(work), '<WORK>'), (str(SET/'inputs'), '<INPUTS>'),
             (str(TOOLS), '<TOOLS>'), (str(BINARY), '<BINARY>'), (str(REPO), '<REPO>')]
    told = dict(returncode=done.returncode, summaries={})
    for path in sorted(out.rglob('summary.json')):
        told['summaries'][str(path.relative_to(out))] = normalized(json.loads(path.read_text()), swaps)
    return told


def main(argv):
    if not argv or argv[0].startswith('-'):
        print(__doc__); return 2
    work = Path(argv[0]).resolve(); work.mkdir(parents=True, exist_ok=True)
    (work/'cache').mkdir(exist_ok=True)
    only = argv[argv.index('--only')+1] if '--only' in argv else None
    update = '--update' in argv
    if only is None:                       # the cache cases share a file, in order
        for stale in (work/'cache').glob('*'):
            stale.unlink()
    cases = [c for c in json.loads((SET/'cases.json').read_text()) if only in (None, c['name'])]
    failed = 0
    for case in cases:
        start = time.time()
        told = run(case, work)
        golden = SET/'golden'/f"{case['name']}.json"
        if update:
            golden.write_text(json.dumps(told, indent=1, sort_keys=True)+'\n')
            print(f"{case['name']:30} written   {time.time()-start:5.1f}s"); continue
        if not golden.is_file():
            print(f"{case['name']:30} NO GOLDEN FILE"); failed += 1; continue
        found = list(differences(timed_out_alike(json.loads(golden.read_text())), timed_out_alike(json.loads(json.dumps(told)))))
        print(f"{case['name']:30} {'ok' if not found else 'DIFFERS':8} {time.time()-start:5.1f}s")
        for line in found[:8]:
            print('    '+line)
        failed += bool(found)
    print(f'{failed} of {len(cases)} cases differ' if not update else f'{len(cases)} golden files written')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
