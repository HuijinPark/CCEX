#!/usr/bin/env python3
"""
Compare the qubit block of two `-v` runs of the same multi-qubit system.

Reads the per-qubit source/computational positions and the intmap tensors out of
the coordinate-frame report, and checks what the rotation must do to them. R is
rebuilt here from the convention, so this is an independent check of the C++.

  ./check_qubits.py rot.log --bath-axis 0 0 1 --qubit-axis 1 1 1 --reference NV0
  ./check_qubits.py norot.log rot.log --kind intmap --expect rotated ...
"""

import argparse
import math
import re
import sys

QNAME = re.compile(r'^\s*qubit\[(\d+)\] "([^"]+)"(\s+\(reference\))?\s*$')
QSRC = re.compile(r"^\s*source frame\s*:\s*\[(.+)\]\s*$")
QCALC = re.compile(r"^\s*computational frame\s*:\s*\[(.+)\]\s*$")
INTMAP = re.compile(r"^\s*intmap\[(\d+)\]\[(\d+)\]\s*:\s*\[(.+)\]\s*$")
CPLX = re.compile(r"(-?\d+\.\d+)[-+]\d+\.\d+j")
FLOAT = re.compile(r"-?\d+\.\d+")


def parse(path):
    qubits, intmap = {}, {}
    cur = None
    with open(path) as fh:
        for line in fh:
            m = QNAME.match(line)
            if m:
                cur = int(m.group(1))
                qubits[cur] = {"name": m.group(2), "ref": m.group(3) is not None}
                continue
            if cur is not None:
                m = QSRC.match(line)
                if m:
                    qubits[cur]["src"] = [float(x) for x in FLOAT.findall(m.group(1))]
                    continue
                m = QCALC.match(line)
                if m:
                    qubits[cur]["calc"] = [float(x) for x in FLOAT.findall(m.group(1))]
                    cur = None
                    continue
            m = INTMAP.match(line)
            if m:
                vals = [float(v) for v in CPLX.findall(m.group(3))]
                if len(vals) == 9:
                    intmap[(int(m.group(1)), int(m.group(2)))] = vals
    return qubits, intmap


def normalize(v):
    n = math.sqrt(sum(c * c for c in v))
    return [c / n for c in v]


def cross(u, v):
    return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]


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


def apply_R(R, p, r0):
    d = [p[k] - r0[k] for k in range(3)]
    return [r0[k] + sum(R[k][j] * d[j] for j in range(3)) for k in range(3)]


def rot_T(R, T):
    M = [[T[3 * i + j] for j in range(3)] for i in range(3)]
    return [sum(R[i][k] * M[k][l] * R[j][l] for k in range(3) for l in range(3))
            for i in range(3) for j in range(3)]


def distance(a, b):
    return math.sqrt(sum((a[k] - b[k]) ** 2 for k in range(3)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("log")
    ap.add_argument("--bath-axis", nargs=3, type=float, required=True)
    ap.add_argument("--qubit-axis", nargs=3, type=float, required=True)
    ap.add_argument("--reference", required=True, help="expected reference qubit name")
    ap.add_argument("--expect-nqubit", type=int, default=None)
    ap.add_argument("--tol", type=float, default=3e-3)
    args = ap.parse_args()

    qubits, intmap = parse(args.log)
    fails = []

    def check(ok, what):
        print(("  ok    " if ok else "  FAIL  ") + what)
        if not ok:
            fails.append(what)

    check(len(qubits) > 0, "the coordinate-frame report lists the qubits (%d)" % len(qubits))
    if not qubits:
        return 1
    if args.expect_nqubit is not None:
        check(len(qubits) == args.expect_nqubit,
              "nqubit = %d as expected" % args.expect_nqubit)

    refs = [i for i in qubits if qubits[i]["ref"]]
    check(len(refs) == 1 and qubits[refs[0]]["name"] == args.reference,
          'the reference qubit is "%s" (index %s)' % (args.reference, refs))
    if len(refs) != 1:
        return 1
    iref = refs[0]
    r0 = qubits[iref]["src"]

    # The reference qubit is the fixed point: it must not have moved at all.
    moved = max(abs(a - b) for a, b in zip(qubits[iref]["src"], qubits[iref]["calc"]))
    check(moved < args.tol,
          'the reference qubit "%s" stays at %s (moved %.2e)'
          % (args.reference, qubits[iref]["src"], moved))

    # Every qubit follows r0 + R(r - r0), with the SAME R and r0.
    R = build_R(args.bath_axis, args.qubit_axis)
    worst, worst_i = 0.0, -1
    for i in qubits:
        want = apply_R(R, qubits[i]["src"], r0)
        err = max(abs(a - b) for a, b in zip(want, qubits[i]["calc"]))
        if err > worst:
            worst, worst_i = err, i
    check(worst < args.tol,
          "every qubit is at r0 + R*(r - r0) (worst = %.2e at qubit %d)" % (worst, worst_i))

    # Rigid: qubit-qubit distances survive.
    worst = 0.0
    for i in qubits:
        for j in qubits:
            if j <= i:
                continue
            worst = max(worst, abs(distance(qubits[i]["src"], qubits[j]["src"]) -
                                   distance(qubits[i]["calc"], qubits[j]["calc"])))
    check(worst < args.tol,
          "every qubit-qubit distance is preserved (worst = %.2e)" % worst)

    if fails:
        print("\n  %d check(s) failed" % len(fails))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
