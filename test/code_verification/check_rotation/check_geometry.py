#!/usr/bin/env python3
"""
Compare an unrotated CCEX run against a rotated one and check what the rotation
is supposed to leave alone.

Both logs must come from `-v` runs of the SAME system, differing only in the
bath_coordinate_rotation section. Positions are read out of the BathArray report,
which prints 3 decimals -- that sets the tolerances below, nothing else.

R is rebuilt here from the convention (ez = q, ex = q x b, ey = ez x ex, as ROWS)
rather than being hard-coded, so this is an independent check of the C++ and not a
copy of its output.

  ./check_geometry.py norot.log rot.log --bath-axis 0 0 1 --qubit-axis 1 1 1 --r0 10 20 30
"""

import argparse
import math
import re
import sys

# "      [  0]   13C  12.000  20.000  30.000 ( S = 0.5, gyro =      6.728, mainspidx = 0, r_qb =   2.000 )"
BATHLINE = re.compile(
    r"^\s*\[\s*(\d+)\]\s+(\S+)\s+"
    r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+"
    r"\(\s*S\s*=\s*(-?\d+\.\d+),\s*gyro\s*=\s*(-?\d+\.\d+),\s*"
    r"mainspidx\s*=\s*(-?\d+),\s*r_qb\s*=\s*(-?\d+\.\d+)\s*\)")

# "      xyz (A)           :   [ 10.00     , 20.00     , 30.00       ]"
QUBITXYZ = re.compile(
    r"^\s*xyz \(A\)\s*:\s*\[\s*(-?\d+\.\d+)\s*,\s*(-?\d+\.\d+)\s*,\s*(-?\d+\.\d+)\s*\]")

POS_TOL  = 3e-3   # 3 printed decimals
DIST_TOL = 5e-3


def parse(path):
    """-> (list of spin records ordered by index, qubit xyz)

    Both the per-spin line and the qubit xyz are printed more than once -- the spins
    while each bath file is read and again in the final BathArray report, the qubit
    once before the files are read (still [0,0,0]) and again after. The LAST
    occurrence of each is the final state, so later lines overwrite earlier ones."""
    spins, qubit = {}, None
    with open(path) as fh:
        for line in fh:
            m = BATHLINE.match(line)
            if m:
                spins[int(m.group(1))] = {
                    "name": m.group(2),
                    "xyz": (float(m.group(3)), float(m.group(4)), float(m.group(5))),
                    "spin": float(m.group(6)),
                    "gyro": float(m.group(7)),
                    "mainspidx": int(m.group(8)),
                    "r_qb": float(m.group(9)),
                }
                continue
            m = QUBITXYZ.match(line)
            if m:
                qubit = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
    return [spins[i] for i in sorted(spins)], qubit


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
    return [ex, ey, ez]          # rows


def apply_R(R, p, r0):
    d = [p[k] - r0[k] for k in range(3)]
    return [r0[k] + sum(R[k][j] * d[j] for j in range(3)) for k in range(3)]


def distance(a, b):
    return math.sqrt(sum((a[k] - b[k]) ** 2 for k in range(3)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("before")
    ap.add_argument("after")
    ap.add_argument("--bath-axis", nargs=3, type=float, required=True)
    ap.add_argument("--qubit-axis", nargs=3, type=float, required=True)
    ap.add_argument("--r0", nargs=3, type=float, required=True)
    ap.add_argument("--expect-nspin", type=int, default=None)
    # The qubit check below reads the LAST "xyz (A)" in the log, which is the last
    # qubit -- fine for one qubit, wrong for several. check_qubits.py covers every
    # qubit properly, so multi-qubit runs turn this one off rather than mis-read it.
    ap.add_argument("--skip-qubit-check", action="store_true")
    args = ap.parse_args()

    a, qa = parse(args.before)
    b, qb = parse(args.after)

    fails = []

    def check(ok, what):
        print(("  ok    " if ok else "  FAIL  ") + what)
        if not ok:
            fails.append(what)

    check(len(a) > 0 and len(a) == len(b),
          "same number of bath spins survive the rbath cut (%d vs %d)" % (len(a), len(b)))
    if len(a) != len(b) or not a:
        print("\n  cannot continue: spin lists do not line up")
        return 1

    if args.expect_nspin is not None:
        check(len(a) == args.expect_nspin,
              "nspin = %d as expected (all bath files were read)" % args.expect_nspin)

    # [9] order, labels and per-spin bookkeeping
    check([s["name"] for s in a] == [s["name"] for s in b],
          "bath-spin labels are unchanged and in the same order")
    check([s["spin"] for s in a] == [s["spin"] for s in b] and
          [s["gyro"] for s in a] == [s["gyro"] for s in b],
          "bath-spin species properties (spin, gyro) are unchanged")
    check([s["mainspidx"] for s in a] == [s["mainspidx"] for s in b],
          "mainspidx assignments are unchanged")

    # [6] the qubit does not move
    r0 = args.r0
    if args.skip_qubit_check:
        print("  --    qubit position check skipped (multi-qubit; see check_qubits.py)")
    else:
        check(qa is not None and qb is not None and
              max(abs(qb[k] - r0[k]) for k in range(3)) < 1e-2 and
              max(abs(qa[k] - qb[k]) for k in range(3)) < 1e-2,
              "the central qubit stays at [%g, %g, %g] in both runs" % tuple(r0))

    # [7] every qubit-bath distance
    worst = max(abs(x["r_qb"] - y["r_qb"]) for x, y in zip(a, b))
    check(worst < DIST_TOL,
          "every qubit-bath distance is preserved (worst = %.2e)" % worst)

    # [8] every bath-bath distance
    worst = 0.0
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            worst = max(worst, abs(distance(a[i]["xyz"], a[j]["xyz"]) -
                                   distance(b[i]["xyz"], b[j]["xyz"])))
    check(worst < DIST_TOL,
          "every bath-bath distance is preserved (worst = %.2e, %d pairs)"
          % (worst, len(a) * (len(a) - 1) // 2))

    # the positions themselves : r_new = r0 + R*(r_old - r0), R rebuilt here
    R = build_R(args.bath_axis, args.qubit_axis)
    worst, worst_i = 0.0, -1
    for i, (x, y) in enumerate(zip(a, b)):
        exp = apply_R(R, x["xyz"], r0)
        err = max(abs(exp[k] - y["xyz"][k]) for k in range(3))
        if err > worst:
            worst, worst_i = err, i
    check(worst < POS_TOL,
          "every rotated position matches r0 + R*(r_old - r0) (worst = %.2e at spin %d)"
          % (worst, worst_i))

    # R itself must send the qubit axis to +z
    q = normalize(args.qubit_axis)
    Rq = [sum(R[k][j] * q[j] for j in range(3)) for k in range(3)]
    check(max(abs(Rq[k] - (1.0 if k == 2 else 0.0)) for k in range(3)) < 1e-12,
          "R * normalize(qubit_axis) = [0, 0, 1]")

    if fails:
        print("\n  %d check(s) failed" % len(fails))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
