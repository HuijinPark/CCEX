#!/usr/bin/env python3
"""Fault-injecting stand-in for the MPI launcher (plan 0A, Tasks 12 and 15). Test only.

Used as --mpiexec. CCEX_FAULT='<step>:<action>[;...]' names a step directory by
its tail relative to the run root ('0005-plan' or 'condition-001/0006-run').
Actions: fail (exit 3, no output), sleep=<s> (then run normally), hang-mpi (real
MPI ranks that sleep, so a timeout has to kill and reap them), ignore-term
(ignore SIGTERM and sleep, so the launcher needs SIGKILL), perturb-curve (run,
then scale the last curve value by 1+1e-3), drop-probe (run, then delete the
rank-0 probe TSVs; an except step writes its excepted curves as the same files),
truncate-probe (run, then keep the header and first row of each rank-0 probe TSV),
drop-factors (run, then delete the rank-0 TSV of a factors step).

CCEX_PROFILE=flat|near|cube rewrites the '[profile]' timings of every call with
values fixed by cluster count and dimension, so measured costs are reproducible:
flat gives every cluster the same time (a candidate's total still scales with its
cluster count, so the cost route rejects every split), near adds a 1e-6 relative
difference per dimension, and cube follows the D^3 work proxy so the cost route
behaves as it does with real measurements. CCEX_FAULT_LOG names a
JSONL file that records every action taken, so a case can check that the
injection fired at the intended step. Every other call execs CCEX_REAL_MPIEXEC
(default mpirun) with the same arguments.
"""
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

FIELD = re.compile(r'(\w+)=(\S+)')
TIMING = re.compile(r'\b(build_s|setup_s|evolution_s)=\S+')
args = sys.argv[1:]
step = Path(args[args.index('-f')+1]).resolve().parent
real = os.environ.get('CCEX_REAL_MPIEXEC', 'mpirun')
mode = os.environ.get('CCEX_PROFILE')
rules = dict(rule.rsplit(':', 1) for rule in os.environ.get('CCEX_FAULT', '').split(';') if rule)
rule, action = next(((s, a) for s, a in rules.items() if step.as_posix().endswith('/'+s)), (None, None))


def record(**fields):
    path = os.environ.get('CCEX_FAULT_LOG')
    if path:
        with open(path, 'a') as log:
            log.write(json.dumps(dict(rule=rule, action=action, step=str(step), **fields))+'\n')


UNIT = dict(flat=lambda d: 1e-3,                 # every cluster costs the same
            near=lambda d: 1e-3*(1+1e-6*d),      # almost the same, ordered by dimension
            cube=lambda d: 1e-9*d**3)            # follows the D^3 work proxy


def rewrite(line):
    """Replace measured timings with values fixed by cluster count and dimension."""
    fields = dict(FIELD.findall(line))
    fixed = dict(build_s=0., setup_s=0., evolution_s=int(fields['clusters'])*UNIT[mode](int(fields['dim'])))
    return TIMING.sub(lambda m: f'{m.group(1)}={fixed[m.group(1)]!r}', line)


def launch():
    if mode is None:
        return subprocess.run([real, *args]).returncode
    done = subprocess.run([real, *args], capture_output=True, text=True)
    sys.stdout.write(''.join((rewrite(line) if line.startswith('[profile] ') else line)+'\n'
                             for line in done.stdout.splitlines()))
    sys.stderr.write(done.stderr)
    return done.returncode


if action is None:
    if mode is None:
        os.execvp(real, [real, *args])
    sys.exit(launch())
if action == 'fail':
    record(); sys.exit(3)
if action == 'ignore-term':
    record(); signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(600); sys.exit(4)
if action == 'hang-mpi':
    marker = f'ccex-hang-{os.getpid()}'
    started = step/'hang-started'
    record(marker=marker, started=str(started))
    # Each rank records its own start, so a timeout case can tell "the ranks ran
    # and were reaped" from "the limit hit before they started".
    program = ("import os, sys, time\n"
               "open('{}.{}'.format(sys.argv[1], os.getpid()), 'w').write('started')\n"
               "time.sleep(600)")
    binary = args.index('-f')-1
    os.execvp(real, [real, *args[:binary], sys.executable, '-c', program, str(started), marker])
record()
if action.startswith('sleep='):
    time.sleep(float(action[len('sleep='):]))
    sys.exit(launch())
code = launch()
if action == 'perturb-curve':
    for f in sorted(step.glob('response.curve.state-*.tsv')):
        rows = f.read_text().splitlines()
        probe, t, real_part, imag = rows[-1].split('\t')
        rows[-1] = '\t'.join((probe, t, repr(float(real_part)*(1+1e-3)), imag))
        f.write_text('\n'.join(rows)+'\n')
elif action == 'drop-probe':
    for f in step.glob('response.rank-0.state-*.tsv'):
        f.unlink()
elif action == 'drop-factors':
    for f in step.glob('response.rank-0.factors.tsv'):
        f.unlink()
elif action == 'truncate-probe':
    for f in step.glob('response.rank-0.state-*.tsv'):
        f.write_text('\n'.join(f.read_text().splitlines()[:2])+'\n')
else:
    sys.exit(f'unknown fault action {action}')
sys.exit(code)
