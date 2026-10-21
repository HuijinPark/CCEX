#!/usr/bin/env python3
"""Method D in C++ (--regroup) against its prototype, on the stored P1 probes ((internal notes),
stages 0 and 4). Needs the stored probes of ~/tmp/ccex-p1-cluster, so it is a study, not a standing check.

    rebundle2.py trace 51-boundary-pcce.json > PY.jsonl        ((internal notes); with
                                                              REBUNDLE_NOREVISIT=1, and REBUNDLE_SWAP=1 for the swapped split)
    cpp/bin/regroup_replay PY.jsonl > CPP.jsonl                  (the C++ decisions on the same probes)
    regroup_parity.py PY.jsonl CPP.jsonl                         (the comparison; exit 1 when a case differs)

A case is the same when every decision of every round is: the partition, the worst point (s*, t*), U, the culprit pairs in
order, each region (skipped or not, its neighbours, its candidate count, the rejections by reason, the accepted candidate),
the accepted partition, and at the end the stop reason, the final partition, the changes and the regions examined. Numbers
(excess, E_loc of the region, of each neighbour before and after, of the accepted candidate) may differ by 1e-9; a curve may
differ by 1e-9 (complex absolute difference) at the points finite on both sides and has to be non-finite at the same points.
The first difference of a case is reported with its path."""
import json
import math
import sys

TOL = 1e-9
EXACT = ('units', 's', 't', 'U', 'culprits', 'accepted')
REGION_EXACT = ('culprit', 'region', 'skipped', 'neighbours', 'candidates', 'rejected', 'accepted', 'U_new')
REGION_NUMBERS = ('e_now', 'e', 'excess_new')
REGION_LISTS = ('before', 'after')


def number(a, b):
    if isinstance(a, str) or isinstance(b, str):
        return a == b
    return abs(a - b) <= TOL


def curve(a, b, at):
    for s, (x, y) in enumerate(zip(a, b)):
        for t, (p, q) in enumerate(zip(x, y)):
            fp = all(not isinstance(v, str) for v in p)
            fq = all(not isinstance(v, str) for v in q)
            if fp != fq:
                return f'{at}[{s}][{t}]: non-finite on one side only ({p} vs {q})'
            if fp and math.hypot(p[0] - q[0], p[1] - q[1]) > TOL:
                return f'{at}[{s}][{t}]: |diff| {math.hypot(p[0] - q[0], p[1] - q[1]):.3g}'
    if len(a) != len(b) or any(len(x) != len(y) for x, y in zip(a, b)):
        return f'{at}: shapes differ'
    return None


def compare(py, cpp):
    """the first difference between two traces, or None"""
    if len(py) != len(cpp):
        return f'rounds: {len(py) - 1} in the prototype, {len(cpp) - 1} in C++'
    for k, (a, b) in enumerate(zip(py, cpp)):
        if 'stop' in a:
            for key in ('stop', 'final', 'changes', 'examined'):
                if a[key] != b.get(key):
                    return f'end/{key}: {a[key]!r} vs {b.get(key)!r}'
            continue
        for key in EXACT:
            if a[key] != b.get(key):
                return f'round {k}/{key}: {a[key]!r} vs {b.get(key)!r}'
        if not number(a['excess'], b['excess']):
            return f'round {k}/excess: {a["excess"]!r} vs {b["excess"]!r}'
        found = curve(a['curve'], b['curve'], f'round {k}/curve')
        if found:
            return found
        if len(a['regions']) != len(b['regions']):
            return f'round {k}/regions: {len(a["regions"])} vs {len(b["regions"])}'
        for r, (x, y) in enumerate(zip(a['regions'], b['regions'])):
            at = f'round {k}/region {r}'
            for key in REGION_EXACT:
                if x.get(key) != y.get(key):
                    return f'{at}/{key}: {x.get(key)!r} vs {y.get(key)!r}'
            for key in REGION_NUMBERS:
                if (key in x) != (key in y) or (key in x and not number(x[key], y[key])):
                    return f'{at}/{key}: {x.get(key)!r} vs {y.get(key)!r}'
            for key in REGION_LISTS:
                if (key in x) != (key in y) or (key in x and (len(x[key]) != len(y[key]) or
                                                             not all(map(number, x[key], y[key])))):
                    return f'{at}/{key}: {x.get(key)!r} vs {y.get(key)!r}'
    return None


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    key = lambda d: (d['CR'], d['field'], d['cfg'], d['N'], d['K'])
    py = {key(d): d for d in map(json.loads, open(argv[0]))}
    cpp = {key(d): d for d in map(json.loads, open(argv[1]))}
    if set(py) != set(cpp):
        print(f'cases differ: {len(py)} in the prototype, {len(cpp)} in C++')
        return 1
    same, differ, worst = 0, [], 0.0
    for k in sorted(py):
        found = compare(py[k]['trace'], cpp[k]['trace'])
        if found:
            differ.append((k, found))
        else:
            same += 1
        for a, b in zip(py[k]['trace'], cpp[k]['trace']):
            if 'curve' in a and 'curve' in b:
                for x, y in zip(a['curve'], b['curve']):
                    for p, q in zip(x, y):
                        if not any(isinstance(v, str) for v in p + q):
                            worst = max(worst, math.hypot(p[0] - q[0], p[1] - q[1]))
    changed = sum(1 for k in py if py[k]['trace'][-1]['changes'])
    print(f'{len(py)} cases: {same} identical, {len(differ)} different; {changed} with a change in the prototype; '
          f'largest curve difference at a finite point {worst:.3g} (where the rounds match)')
    for k, found in differ:
        print(f'  {"/".join(map(str, k))}: {found}')
    return 1 if differ else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
