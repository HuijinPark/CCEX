#!/usr/bin/env python3
"""compare.py RUN_DIR REF_DIR [--expected expected.json] [--tol 1e-9]

For every case, the largest |run - reference| over the files both wrote. The verdict is taken
on the _noDiv files and the density-matrix files; the _wiDiv coherence (normalised by the
free evolution) is printed but not judged, because it amplifies rounding wherever the
denominator is small. A case listed in expected.json for this reference is reported as
"expected" with its reason instead of FAIL. A case with no output, no reference, a reference file the
run did not write, nothing judged in common, or an output that is empty, of another length or, where judged,
non-finite where the reference is finite, fails as well, whatever expected.json says; a file the
reference wrote under another name is compared under the name in the reference's "renames"
({case: {run file: reference file}}). No case at all fails. Exit status 1 if any case fails."""
import sys, os, glob, json, numpy as np
args = sys.argv[1:]; run, ref = args[0], args[1]
tol = float(args[args.index('--tol') + 1]) if '--tol' in args else 1e-9
expected = json.load(open(args[args.index('--expected') + 1])) if '--expected' in args else {}
version = os.path.basename(ref.rstrip('/')); exp = expected.get(version, {}); tol = float(exp.get('tol', tol))
renames = exp.get('renames', {})
def load(path):
    """Every whitespace-separated token as a number (complex() reads the whole token, inf and nan included), and
    the first token it cannot read, as "line N: token", or None. Such a file is broken whatever the reference holds
    there. Until 2026-10-02 a pattern read each line and took "1.0+infj" as 1.0 (review 78d583 R2), then for a while
    rescued "1.0BROKEN" as 1.0 (7eed38 R2), then read it as a nan, which passed beside a nan of the reference (fa7f39 R1)."""
    out, bad = [], None
    for n, line in enumerate(open(path), 1):
        for tok in line.split():
            try: out.append(complex(tok))
            except ValueError:
                out.append(complex('nan')); bad = bad or f'line {n}: "{tok}"'
    return np.array(out), bad
def broken(a, b, judged):
    """Why this output cannot stand against its reference whatever expected.json says (review 1bd0ff R1):
    empty, a different length, or, in a judged file, a non-finite value where the reference is finite. A
    non-finite value only in the reference (the nan of the old dm output) is a numerical difference."""
    if a.size == 0: return 'empty'
    if a.shape != b.shape: return f'{a.size} values against {b.size}'
    if judged and (~np.isfinite(a) & np.isfinite(b)).any(): return 'a non-finite value the reference does not have'
    return None
def worst(a, b):
    if a.shape != b.shape: return float('inf')
    d = np.abs(a - b); ok = np.isfinite(d)
    if (~ok & ~(np.isnan(a) & np.isnan(b))).any(): return float('inf')      # a nan on one side only
    return float(d[ok].max()) if ok.any() else 0.0
fails = 0; rows = []
# The cases are those of the run and of the reference together, and a case fails when nothing could
# be compared: no output, no reference, a reference file the run did not write, or no judged file in
# common. Until 2026-10-02 these were notes, and a run that wrote nothing passed (review R1).
names = sorted({os.path.basename(d.rstrip('/')) for d in glob.glob(f'{run}/*/') + glob.glob(f'{ref}/*/')})
for name in names:
    cdir, rdir = f'{run}/{name}/', f'{ref}/{name}'
    files = sorted(f for f in os.listdir(cdir) if f != 'run.log') if os.path.isdir(cdir) else []
    if not files or not os.path.isdir(rdir): rows.append((name, None, None, 'no reference' if files else 'no output')); continue
    judged, shown, notes, compared, missing = 0.0, 0.0, [], 0, []
    names_in_ref = {f: renames.get(name, {}).get(f, f) for f in files}    # run file -> reference file
    for f in files:
        rf = names_in_ref[f]
        if not os.path.exists(f'{rdir}/{rf}'): notes.append(f'only here: {f}'); continue
        if rf != f: notes.append(f'{f} compared with {rf}')
        (a, bad_a), (b, bad_b) = load(cdir + f), load(f'{rdir}/{rf}')
        why = (f'unreadable {bad_a}' if bad_a else f'reference unreadable {bad_b}' if bad_b else
               broken(a, b, not f.endswith('_wiDiv')))
        if why: notes.insert(0, f'broken {f}: {why}'); missing.append(f); continue
        w = worst(a, b)
        if f.endswith('_wiDiv'): shown = max(shown, w)
        else: judged = max(judged, w); compared += 1
    for f in sorted(os.listdir(rdir)):
        if f != 'run.log' and f not in names_in_ref.values(): notes.append(f'only in reference: {f}'); missing.append(f)
    if missing or not compared: rows.append((name, None, shown, '; '.join(['not compared in full'] + notes))); continue    # also a broken file
    rows.append((name, judged, shown, '; '.join(notes)))
print(f'{"case":24} {"noDiv/dm":>10} {"wiDiv":>10}  verdict   (reference {version}, tolerance {tol:g})')
for name, j, s, note in rows:
    reason = exp.get('cases', {}).get(name)
    if j is None:     # nothing or not everything compared: no reason excuses it (reviews 4cc510 R2, 8fecd5 R1)
        fails += 1; print(f'{name:24} {"-":>10} {"-":>10}  {"FAIL":8}  {note}'); continue
    if j <= tol: verdict = 'ok'
    elif reason: verdict = 'expected'
    else: verdict = 'FAIL'; fails += 1
    print(f'{name:24} {j:10.3e} {s:10.3e}  {verdict:8}  {reason or note}')
if not names: print(f'no case in {run} or {ref}'); fails += 1
print(f'{fails} case(s) beyond {tol:g} without an explanation')
sys.exit(1 if fails else 0)
