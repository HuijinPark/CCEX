#!/usr/bin/env python3
"""
Check that the stored 3x3 tensors moved with the frame: T_new = R * T_old * R^T.

Reads the hypf[i][j] and quad[i] lines out of two `-v` logs of the same system --
one unrotated, one rotated -- and compares them index by index. R is rebuilt here
from the convention rather than read out of the log, so this is an independent
check of the C++.

Matching index by index is also what makes this a test of the tensor-file LOOKUP:
the tensor files are keyed by position and are read before anything moves, so if
the lookup had used rotated coordinates instead, spin i would have picked up some
other atom's tensor and no rotation of the old value could reproduce it.

CCEX prints tensor components with "%3.2f", so the comparison is absolute and the
tolerance is set by that, not by the arithmetic.

  ./check_tensors.py norot.log rot.log --bath-axis 0 0 1 --qubit-axis 1 1 1 \\
                     --expect rotated --kind hypf
"""

import argparse
import math
import re
import sys

HYPF = re.compile(r"^\s*hypf\[(\d+)\]\[(\d+)\]\s*:\s*\[(.+)\]\s*$")
QUAD = re.compile(r"^\s*quad\[(\d+)\]\s*:\s*\[(.+)\]\s*$")
INTMAP = re.compile(r"^\s*intmap\[(\d+)\]\[(\d+)\]\s*:\s*\[(.+)\]\s*$")
# "1.46-0.00j" / "-2.04+0.00j"
CPLX = re.compile(r"(-?\d+\.\d+)[-+]\d+\.\d+j")


def parse(path, kind):
    """-> {key: [9 real components]}; later occurrences win."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if kind == "hypf":
                m = HYPF.match(line)
            elif kind == "intmap":
                m = INTMAP.match(line)
            else:
                m = QUAD.match(line)
            if m is None:
                continue
            if kind in ("hypf", "intmap"):
                key, body = (int(m.group(1)), int(m.group(2))), m.group(3)
            else:
                key, body = int(m.group(1)), m.group(2)
            vals = [float(v) for v in CPLX.findall(body)]
            if len(vals) == 9:
                out[key] = vals
    return out


def normalize(v):
    n = math.sqrt(sum(c * c for c in v))
    return [c / n for c in v]


def cross(u, v):
    return [u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0]]


def build_R(bath_axis, qubit_axis):
    b, q = normalize(bath_axis), normalize(qubit_axis)
    c = cross(q, b)
    # Collinear axes leave the azimuth of the new x undefined. The C++ accepts exactly
    # one such case -- both axes on +z, which is the identity -- and refuses the rest;
    # mirror that here so the checkers agree with what they are checking.
    if math.sqrt(sum(x * x for x in c)) < 1e-12:
        if abs(q[0]) < 1e-12 and abs(q[1]) < 1e-12 and q[2] > 0 and \
           abs(b[0]) < 1e-12 and abs(b[1]) < 1e-12 and b[2] > 0:
            return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        raise SystemExit("bath_axis and qubit_axis are collinear; R is undefined")
    ez = q
    ex = normalize(c)
    ey = cross(ez, ex)
    return [ex, ey, ez]


def rotate(R, T):
    """R * T * R^T, with T flattened row-major."""
    M = [[T[3 * i + j] for j in range(3)] for i in range(3)]
    out = []
    for i in range(3):
        for j in range(3):
            out.append(sum(R[i][k] * M[k][l] * R[j][l]
                           for k in range(3) for l in range(3)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("before")
    ap.add_argument("after")
    ap.add_argument("--kind", choices=["hypf", "quad", "intmap"], required=True)
    ap.add_argument("--bath-axis", nargs=3, type=float, required=True)
    ap.add_argument("--qubit-axis", nargs=3, type=float, required=True)
    ap.add_argument("--expect", choices=["rotated", "unchanged", "mixed"], default="rotated")
    ap.add_argument("--tol", type=float, default=0.05)
    args = ap.parse_args()

    a = parse(args.before, args.kind)
    b = parse(args.after, args.kind)
    fails = []

    def check(ok, what):
        print(("  ok    " if ok else "  FAIL  ") + what)
        if not ok:
            fails.append(what)

    check(len(a) > 0, "%s tensors found in the unrotated log (%d)" % (args.kind, len(a)))
    check(set(a) == set(b),
          "the same %d %s entries appear in both runs" % (len(a), args.kind))
    if not a or set(a) != set(b):
        return 1

    R = build_R(args.bath_axis, args.qubit_axis)

    worst_rot, worst_key = 0.0, None      # |T_after - R T_before R^T|
    worst_same, moved = 0.0, 0            # |T_after - T_before|
    for k in sorted(a, key=str):
        want = rotate(R, a[k]) if args.expect == "rotated" else a[k]
        err = max(abs(x - y) for x, y in zip(want, b[k]))
        if err > worst_rot:
            worst_rot, worst_key = err, k
        same = max(abs(x - y) for x, y in zip(a[k], b[k]))
        worst_same = max(worst_same, same)
        if same > args.tol:
            moved += 1

    if args.expect == "mixed":
        # hf_tensor_frame = "qubit" with readmode 2/3: what the tensor file supplied is
        # already in the computational basis and must be left alone, while the
        # point-dipole fallback for a spin the file did not cover was built from the
        # source geometry and must be transformed exactly once. Both kinds live in the
        # same BathArray, so every entry has to match ONE of the two -- and both classes
        # have to be non-empty, or the test is not exercising the split it claims to.
        kept, rotated, bad = [], [], []
        for k in sorted(a, key=str):
            e_keep = max(abs(x - y) for x, y in zip(a[k], b[k]))
            e_rot = max(abs(x - y) for x, y in zip(rotate(R, a[k]), b[k]))
            if e_keep < args.tol and e_rot < args.tol:
                kept.append(k)          # indistinguishable (e.g. an isotropic tensor)
            elif e_keep < args.tol:
                kept.append(k)
            elif e_rot < args.tol:
                rotated.append(k)
            else:
                bad.append((k, e_keep, e_rot))

        check(not bad,
              "every %s tensor is either left as read or transformed exactly once "
              "(%d unexplained%s)"
              % (args.kind, len(bad),
                 "" if not bad else ", worst %s keep=%.3g rot=%.3g" % bad[0]))
        check(len(kept) > 0,
              "%d %s tensors came from the file and were left alone" % (len(kept), args.kind))
        check(len(rotated) > 0,
              "%d %s tensors were point-dipole fallbacks and were transformed"
              % (len(rotated), args.kind))
        print("        matched (kept) : %s ..." % [str(k) for k in kept[:4]])
        print("        fallback (rot) : %s ..." % [str(k) for k in rotated[:4]])

    elif args.expect == "rotated":
        check(worst_rot < args.tol,
              "every %s tensor satisfies T_new = R*T_old*R^T (worst = %.3g at %s, tol %.3g)"
              % (args.kind, worst_rot, worst_key, args.tol))
        # If nothing actually changed, the check above would pass for R = I and prove
        # nothing about the transform.
        check(moved > 0,
              "%d of %d %s tensors genuinely changed (the transform is not a no-op)"
              % (moved, len(a), args.kind))
    else:
        check(worst_rot < args.tol,
              "every %s tensor is left exactly as read (worst = %.3g, tol %.3g)"
              % (args.kind, worst_rot, args.tol))

    # A similarity transform by an orthogonal matrix preserves these, whichever branch
    # was taken -- a useful cross-check that the parse itself lined up. It holds in the
    # mixed case too, since both branches are similarity transforms (one of them by I).
    worst_tr, worst_fro = 0.0, 0.0
    for k in a:
        tr_a = a[k][0] + a[k][4] + a[k][8]
        tr_b = b[k][0] + b[k][4] + b[k][8]
        worst_tr = max(worst_tr, abs(tr_a - tr_b))
        fa = math.sqrt(sum(x * x for x in a[k]))
        fb = math.sqrt(sum(x * x for x in b[k]))
        worst_fro = max(worst_fro, abs(fa - fb))
    check(worst_tr < args.tol and worst_fro < args.tol,
          "trace and Frobenius norm are invariant (worst %.3g / %.3g)" % (worst_tr, worst_fro))

    if fails:
        print("\n  %d check(s) failed" % len(fails))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
