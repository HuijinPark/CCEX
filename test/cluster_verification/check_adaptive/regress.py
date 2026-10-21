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
comparison against Python that this replaces. A case whose bin/main.out runs came from a
binary built with MKL (their run.log says "MKL for large clusters : on"; A46) keeps the
decisions exact and lets a float move by 1e-5 of the golden value (close() below).

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


MKL_ON = 'MKL for large clusters : on'      # the banner line of a binary built with MKL (A46)


def built_with_mkl(out):
    """Whether the bin/main.out runs of this case came from a binary built with MKL: their run logs say so."""
    return any(MKL_ON in p.read_text(errors='replace') for p in out.rglob('*.log'))


def close(a, b, mkl):
    """A float of golden (a) against this run's (b). Without MKL: within 1e-9 relative (to max(1, |a|, |b|)). With MKL
    (A46, user decision 2026-10-13): the library moves the last digits of the large clusters, and the policies' error
    numbers (differences of nearly equal coherences) carry that up to 4.4e-6 relative (WSL, 8 of 33 cases); so within
    1e-5 of |a|, and never tighter than 1e-12 absolute (a = 0, and round-off-sized values: regroup's
    max_abs_difference_from_prediction is 1.0e-15 in golden and 8.9e-16 with MKL). Decisions - every int, bool and
    string - stay exact either way."""
    if a == b or (math.isnan(a) and math.isnan(b)):
        return True
    if mkl:
        return abs(a-b) <= max(1e-5*abs(a), 1e-12)
    return abs(a-b) <= 1e-9*max(1., abs(a), abs(b))


def differences(a, b, at='', mkl=False):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                yield f'{at}/{k}: only in {"golden" if k in a else "this run"}'
            else:
                yield from differences(a[k], b[k], f'{at}/{k}', mkl)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            yield f'{at}: {len(a)} entries in golden, {len(b)} in this run'
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                yield from differences(x, y, f'{at}[{i}]', mkl)
    elif isinstance(a, float) or isinstance(b, float):
        ok = isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool) \
            and not isinstance(b, bool) and close(a, b, mkl)
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
    return told, built_with_mkl(out)


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
        told, mkl = run(case, work)
        golden = SET/'golden'/f"{case['name']}.json"
        if update:
            golden.write_text(json.dumps(told, indent=1, sort_keys=True)+'\n')
            print(f"{case['name']:30} written   {time.time()-start:5.1f}s"); continue
        if not golden.is_file():
            print(f"{case['name']:30} NO GOLDEN FILE"); failed += 1; continue
        found = list(differences(timed_out_alike(json.loads(golden.read_text())), timed_out_alike(json.loads(json.dumps(told))),
                                 mkl=mkl))
        print(f"{case['name']:30} {'ok' if not found else 'DIFFERS':8} {time.time()-start:5.1f}s"
              + ('  (MKL build: floats within 1e-5 relative)' if mkl else ''))
        for line in found[:8]:
            print('    '+line)
        failed += bool(found)
    print(f'{failed} of {len(cases)} cases differ' if not update else f'{len(cases)} golden files written')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
