#!/usr/bin/env python3
"""The merge policy run three times with the same request - ccex-adaptive in its own
process, main.out --adaptive in ours, and main.out -f with the request written as
an "adaptive" block in the input file - and the three output trees compared.

    inprocess_compare.py WORKDIR

Two things are normalized, both of them the harness's own doing:
  * the output directory's own path, which differs by construction;
  * measured times - the elapsed fields and recorded seconds compare_summaries.py
    drops, and the per-cluster seconds in a [profile] line. Which clusters were
    measured, at which dimension, is compared.
A step's run.log holds CCEX's banner and option report only when the step is
its own process, so only the bracketed lines of a run.log are compared.

Until 2026-10-10 five more cases were here, one per removed policy, and a
single-sample copy of cost-two-conditions (its conditions differed in npulse);
they went with error-gated, budgeted, auto-partition and cost-preserving. The
merge cases after merge-grow were added on 2026-10-11: gCCE with two sampled
states, order 1, prefer-physical and hub candidates. merge-hub-capped is not here:
one of its plan steps is refused (a bundle over bundle.max_dim), which the separate
driver takes as a rejected candidate but which ends the whole process in-process
(MPI_Abort; the same before the 2026-10-11 cleanup).

The regroup cases (method D, --regroup, added 2026-10-11) are here: regroup-p1 changes
the partition once and ends physical, regroup-p1-sparse works on a graph that rdip
leaves sparse, and regroup-p1-limit has three candidates over the expansion limit.
D refuses those before it asks CCEX for anything (adaptive-regroup-plan.md, stage 6),
so the search goes on to the rest of the candidates in-process as it does in the
separate driver, instead of the MPI_Abort a refused plan step ends a process with.
regroup-p1-factors and regroup-p1-except (2026-10-11) take regroup-p1 through the
factors and except work modes (adaptive-regroup-stage3.md, -stage3b.md): steps that
run inside the process too, after the driver's own checks of what CCEX would refuse.
"""
import os
import json, re, shutil, subprocess, sys, time
from pathlib import Path

TEST = Path(__file__).resolve().parent
REPO = TEST.parent.parent.parent
BINARY = Path(os.environ.get('CCEX_BIN', REPO/'bin/main.out'))   # the binary to test, as in regress.py
WORK = Path(sys.argv[1]).resolve(); shutil.rmtree(WORK, ignore_errors=True); WORK.mkdir(parents=True)
CASES = ['merge-grow', 'merge-gcce-s2', 'merge-order1-small', 'merge-physical', 'merge-hub',
         'regroup-p1', 'regroup-p1-sparse', 'regroup-p1-limit', 'regroup-p1-factors', 'regroup-p1-except']
MEASURED = set()      # nothing a policy decides on is priced in measured seconds any more
SECONDS = re.compile(r'(\w+_s)=[^\s]+')     # measured seconds of the [profile] and [phase] lines (the [phase] timers of 864a28774)


def filled(text):
    """regress.py's placeholders; --output and --binary are set again for each way."""
    for key, value in (('<INPUTS>', TEST/'regression/inputs'), ('<WORK>', WORK), ('<TOOLS>', TEST),
                       ('<BINARY>', BINARY), ('<REPO>', REPO)):
        text = text.replace(key, str(value))
    return text


recorded = {c['name']: dict(c, args=[filled(a) for a in c['args']])
            for c in json.loads((TEST/'regression/cases.json').read_text())}


FIXED = ('--input', '--binary', '--cwd', '--output', '--mpiexec', '--ranks')


def as_block(args, out, name):
    """The same request as an "adaptive" object inside a copy of the input file."""
    block, listed = {}, list(args)
    i = 0
    while i < len(listed):
        key = listed[i]
        flag = i+1 == len(listed) or listed[i+1].startswith('--')
        value = True if flag else listed[i+1]
        i += 1 if flag else 2
        if key in FIXED:
            continue
        if not flag:
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    pass
        block[key[2:].replace('-', '_')] = value          # given twice, the last wins, as on the command line
    block['output'] = str(out)
    cfg = json.loads(Path(args[args.index('--input')+1]).read_text())
    cfg['adaptive'] = block
    path = WORK/f'{name}.block.json'
    path.write_text(json.dumps(cfg, indent=2)+'\n')
    return path


# the launcher the other checks use (test/README.md): CCEX_MPI when set (review b9eb09 R2)
MPIRUN = os.environ.get('CCEX_MPI') or 'mpirun'


def run_block(args, out, name):
    start = time.time()
    done = subprocess.run([MPIRUN, '-n', '1', str(BINARY), '-f', str(as_block(args, out, name))],
                          capture_output=True, text=True, cwd=str(REPO))
    return time.time()-start, done.returncode


def run(args, out, inprocess):
    a = list(args)
    a[a.index('--output')+1] = str(out)
    a[a.index('--binary')+1] = str(BINARY)
    # the separate driver launches each step with --mpiexec (default mpirun from PATH): give it the launcher of the
    # in-process run, unless the case names its own (the fault cases) (review 5f7834 R5). The in-process run launches
    # nothing and records its default, so "mpiexec" is the one field the two summaries may differ in for that reason.
    if not inprocess and '--mpiexec' not in a:
        a += ['--mpiexec', MPIRUN]
    cmd = ([MPIRUN, '-n', '1', str(BINARY), '--adaptive'] if inprocess
           else [str(TEST/'cpp/bin/ccex-adaptive')]) + a
    start = time.time()
    done = subprocess.run(cmd, capture_output=True, text=True, cwd=str(REPO))
    return time.time()-start, done.returncode


def normalized(path, tag):
    text = path.read_text(errors='replace').replace(tag, '<OUT>')
    if path.name != 'run.log':
        return text
    return '\n'.join(SECONDS.sub(r'\1=<measured>', line) if line.startswith(('[profile] ', '[phase] ')) else line
                     for line in text.splitlines() if line.startswith('['))


failed = 0
print(f"{'case':34}{'split s':>9}{'in-proc s':>11}{'exit':>8}   status")
for name, args in [(n, recorded[n]['args']) for n in CASES]:
    A, B, C = WORK/f'{name}.split', WORK/f'{name}.inproc', WORK/f'{name}.block'
    ta, ra = run(args, A, False)
    tb, rb = run(args, B, True)
    tc, rc = run_block(args, C, name)
    problems = []
    if ra != rb:
        problems.append(f'exit status {ra} vs {rb}')
    if rc != rb:
        problems.append(f'exit status {rb} on the command line, {rc} from the Adaptive block')
    diff = subprocess.run([sys.executable, str(TEST/'compare_summaries.py'),
                           str(A/'summary.json'), str(B/'summary.json')],
                          capture_output=True, text=True)
    fields = {line.rsplit('/', 1)[-1] for line in diff.stdout.splitlines()}
    launcher = {'mpiexec'} if '--mpiexec' not in args else set()
    if diff.returncode and not fields <= MEASURED | launcher:
        problems.append(f'summary differs in {sorted(fields-MEASURED-launcher)}')
    same = 0
    for one in sorted(A.rglob('*')):
        if not one.is_file() or one.name == 'summary.json':
            continue
        other = B/one.relative_to(A)
        if not other.is_file():
            problems.append(f'missing {one.relative_to(A)}')
        elif normalized(one, f'{name}.split') != normalized(other, f'{name}.inproc'):
            problems.append(f'differs {one.relative_to(A)}')
        else:
            same += 1
    # The block run reads a different input file - the same settings plus the block -
    # so what names or hashes that file differs by construction, and nothing else may.
    own = {'input', 'input_file_sha256', 'input_sha256'}
    diff = subprocess.run([sys.executable, str(TEST/'compare_summaries.py'),
                           str(B/'summary.json'), str(C/'summary.json')], capture_output=True, text=True)
    fields = {line.rsplit('/', 1)[-1] for line in diff.stdout.splitlines()}
    if diff.returncode and not fields <= MEASURED | own:
        problems.append(f'Adaptive block: summary differs in {sorted(fields-MEASURED-own)}')
    for one in sorted(B.rglob('*')):
        if not one.is_file() or one.name == 'summary.json':
            continue
        other = C/one.relative_to(B)
        if not other.is_file():
            problems.append(f'Adaptive block: missing {one.relative_to(B)}')
        elif normalized(one, f'{name}.inproc') != normalized(other, f'{name}.block'):
            problems.append(f'Adaptive block: differs {one.relative_to(B)}')
    status = json.loads((A/'summary.json').read_text()).get('status', '?')
    print(f'{name:34}{ta:9.2f}{tb:11.2f}{ra:>5}/{rb:<3}  {status}  ({same} files alike)')
    for line in problems[:6]:
        print(f'    * {line}')
    failed += bool(problems)
# An input with no "outfile" of its own: the steps name one, and the first of them
# used to copy it into a buffer that only -o or "outfile" allocates (a segfault,
# found on a P1 input on 2026-09-25).
args = recorded['merge-grow']['args']
block = as_block(args, WORK/'no-outfile.block', 'no-outfile')
cfg = json.loads(block.read_text()); cfg.pop('outfile', None); block.write_text(json.dumps(cfg, indent=2)+'\n')
done = subprocess.run([MPIRUN, '-n', '1', str(BINARY), '-f', str(block)],
                      capture_output=True, text=True, cwd=str(REPO))
ok = done.returncode == recorded['merge-grow']['returncode'] and (WORK/'no-outfile.block/summary.json').is_file()
print(f"{'no-outfile (Adaptive block)':34}{'':20}{done.returncode:>5}      {'ok' if ok else 'FAILED'}")
failed += not ok
print(f'{failed} of {len(CASES)+1} cases differ')
sys.exit(1 if failed else 0)
