#!/usr/bin/env python3
"""
Check the magnetic field the run log reports in the source frame.

With the rotation on, "bfield" is read in the COMPUTATIONAL frame, where z is the
qubit axis. R^T takes that back to the source frame, and for a field along the
computational z it has to come out as

    R^T [0,0,Bz] = Bz * normalize(qubit_axis)

which is the sentence "bfield [0,0,500] is 500 G along the physical qubit axis"
written as an equation. The expected value is built here from qubit_axis alone --
R is not used, so this does not just re-derive whatever the code did.

  ./check_srcfield.py run.log --qubit-axis 1 1 1 --bz 500
"""

import argparse
import math
import re
import sys

SRC = re.compile(r"^\s*B field, source frame \(R\^T B\)\s*:\s*\[(.+)\]\s*$")
CALC = re.compile(r"^\s*B field, computational frame\s*:\s*\[(.+)\]\s*$")
FLOAT = re.compile(r"-?\d+\.\d+")


def parse(path):
    src, calc = None, None
    with open(path) as fh:
        for line in fh:
            m = SRC.match(line)
            if m:
                src = [float(x) for x in FLOAT.findall(m.group(1))]
                continue
            m = CALC.match(line)
            if m:
                calc = [float(x) for x in FLOAT.findall(m.group(1))]
    return src, calc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("log")
    ap.add_argument("--qubit-axis", nargs=3, type=float, required=True)
    ap.add_argument("--bz", type=float, required=True)
    ap.add_argument("--tol", type=float, default=1e-6)
    args = ap.parse_args()

    src, calc = parse(args.log)
    fails = []

    def check(ok, what):
        print(("  ok    " if ok else "  FAIL  ") + what)
        if not ok:
            fails.append(what)

    check(calc is not None and src is not None,
          "the report gives the field in both frames")
    if calc is None or src is None:
        return 1

    want_calc = [0.0, 0.0, args.bz]
    check(max(abs(a - b) for a, b in zip(calc, want_calc)) < args.tol,
          "computational-frame field is %s" % (calc,))

    n = math.sqrt(sum(c * c for c in args.qubit_axis))
    want_src = [args.bz * c / n for c in args.qubit_axis]
    worst = max(abs(a - b) for a, b in zip(src, want_src))
    check(worst < args.tol,
          "source-frame field is Bz * normalize(qubit_axis) = %s (worst = %.2e)"
          % ([round(x, 6) for x in want_src], worst))

    # |B| cannot change between frames.
    na = math.sqrt(sum(c * c for c in calc))
    nb = math.sqrt(sum(c * c for c in src))
    check(abs(na - nb) < args.tol,
          "the field magnitude is the same in both frames (%.6f vs %.6f)" % (na, nb))

    if fails:
        print("\n  %d check(s) failed" % len(fails))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
