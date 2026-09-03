#!/usr/bin/env python3
"""
Compare the two sides of the A/B equivalence test.

    A   the same physical system written out by hand in the qubit-aligned frame,
        run with no rotation option
    B   the original input, run with coordinate_frame_rotation on

Both should describe the same physics in the same frame, so everything below is
compared directly -- no transformation is applied here.

Numerical agreement is the criterion, not byte identity: A parses coordinates that
B computes, so the two reach the same numbers by different floating-point routes.
Byte identity of the coherence files is reported when it happens, but never
required.

Tolerances follow what CCEX prints, not what the arithmetic can do:
positions 3 decimals, tensors 2, qubit position 2.

  ./check_ab.py A.log B.log --a-out A --b-out B --kind hypf [--kind quad]
"""

import argparse
import re
import sys

BATHLINE = re.compile(
    r"^\s*\[\s*(\d+)\]\s+(\S+)\s+"
    r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+"
    r"\(\s*S\s*=\s*(-?\d+\.\d+),\s*gyro\s*=\s*(-?\d+\.\d+),\s*"
    r"mainspidx\s*=\s*(-?\d+),\s*r_qb\s*=\s*(-?\d+\.\d+)\s*\)")
HYPF = re.compile(r"^\s*hypf\[(\d+)\]\[(\d+)\]\s*:\s*\[(.+)\]\s*$")
QUAD = re.compile(r"^\s*quad\[(\d+)\]\s*:\s*\[(.+)\]\s*$")
INTMAP = re.compile(r"^\s*intmap\[(\d+)\]\[(\d+)\]\s*:\s*\[(.+)\]\s*$")
# "      xyz (A)           :   [ 10.00     , 20.00     , 30.00       ]"
QUBITXYZ = re.compile(r"^\s*xyz \(A\)\s*:\s*\[(.+)\]\s*$")
CPLX = re.compile(r"(-?\d+\.\d+)[-+]\d+\.\d+j")
BFIELD = re.compile(r"^\s*bfield\s*:\s*\[(.+)\]\s*$")
FLOAT = re.compile(r"-?\d+\.\d+")


def parse(path, kinds):
    spins, tensors, bfield, qxyz = {}, {k: {} for k in kinds}, None, []
    with open(path) as fh:
        for line in fh:
            m = QUBITXYZ.match(line)
            if m:
                v = [float(x) for x in FLOAT.findall(m.group(1))]
                if len(v) == 3:
                    qxyz.append(v)
                continue
            m = BATHLINE.match(line)
            if m:
                spins[int(m.group(1))] = (float(m.group(3)), float(m.group(4)), float(m.group(5)))
                continue
            m = BFIELD.match(line)
            if m:
                v = [float(x) for x in FLOAT.findall(m.group(1))]
                if len(v) == 3:
                    bfield = v
                continue
            if "hypf" in kinds:
                m = HYPF.match(line)
                if m:
                    vals = [float(v) for v in CPLX.findall(m.group(3))]
                    if len(vals) == 9:
                        tensors["hypf"][(int(m.group(1)), int(m.group(2)))] = vals
                    continue
            if "quad" in kinds:
                m = QUAD.match(line)
                if m:
                    vals = [float(v) for v in CPLX.findall(m.group(2))]
                    if len(vals) == 9:
                        tensors["quad"][int(m.group(1))] = vals
                    continue
            if "intmap" in kinds:
                m = INTMAP.match(line)
                if m:
                    vals = [float(v) for v in CPLX.findall(m.group(3))]
                    if len(vals) == 9:
                        tensors["intmap"][(int(m.group(1)), int(m.group(2)))] = vals
    return spins, tensors, bfield, qxyz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a_log")
    ap.add_argument("b_log")
    ap.add_argument("--a-out", required=True, help="A's outfile prefix")
    ap.add_argument("--b-out", required=True, help="B's outfile prefix")
    ap.add_argument("--kind", action="append", default=[], choices=["hypf", "quad", "intmap"])
    # QubitArray_report runs twice, before and after the files are read, so the LAST
    # nqubit entries are the final positions in both runs.
    ap.add_argument("--nqubit", type=int, default=0)
    ap.add_argument("--pos-tol", type=float, default=2e-3)
    ap.add_argument("--tensor-tol", type=float, default=0.02)
    ap.add_argument("--bfield-tol", type=float, default=1e-4)
    ap.add_argument("--coh-tol", type=float, default=1e-8)
    args = ap.parse_args()

    kinds = args.kind or ["hypf"]
    sa, ta, ba, qa_xyz = parse(args.a_log, kinds)
    sb, tb, bb, qb_xyz = parse(args.b_log, kinds)
    fails = []

    def check(ok, what):
        print(("  ok    " if ok else "  FAIL  ") + what)
        if not ok:
            fails.append(what)

    # --- qubit positions ----------------------------------------------------
    if args.nqubit > 0:
        qa_f, qb_f = qa_xyz[-args.nqubit:], qb_xyz[-args.nqubit:]
        ok = len(qa_f) == args.nqubit and len(qb_f) == args.nqubit
        check(ok, "both runs report %d qubit positions" % args.nqubit)
        if ok:
            worst = max(max(abs(x - y) for x, y in zip(p, q)) for p, q in zip(qa_f, qb_f))
            check(worst < 1e-2,   # QubitArray_report prints 2 decimals
                  "qubit positions agree (worst = %.2e)" % worst)

    # --- positions ----------------------------------------------------------
    check(len(sa) > 0 and set(sa) == set(sb),
          "the same %d bath spins, in the same order (A %d / B %d)" % (len(sa), len(sa), len(sb)))
    if not sa or set(sa) != set(sb):
        return 1
    worst = max(max(abs(x - y) for x, y in zip(sa[k], sb[k])) for k in sa)
    check(worst < args.pos_tol,
          "bath positions agree (worst = %.2e, tol %.0e)" % (worst, args.pos_tol))

    # --- tensors ------------------------------------------------------------
    for kind in kinds:
        A, B = ta[kind], tb[kind]
        if not A and not B:
            print("  --    no %s tensors printed in either run" % kind)
            continue
        check(set(A) == set(B), "the same %d %s entries in both runs" % (len(A), kind))
        if set(A) != set(B):
            continue
        worst = max(max(abs(x - y) for x, y in zip(A[k], B[k])) for k in A)
        check(worst < args.tensor_tol,
              "%s tensors agree element by element (worst = %.2e, tol %.0e)"
              % (kind, worst, args.tensor_tol))

    # --- magnetic field ------------------------------------------------------
    if ba is None or bb is None:
        check(False, "Config.bfield was reported in both runs")
    else:
        worst = max(abs(x - y) for x, y in zip(ba, bb))
        check(worst < args.bfield_tol,
              "Config.bfield agrees : A %s vs B %s (worst = %.2e)" % (ba, bb, worst))

    # --- coherence -----------------------------------------------------------
    # The Hamiltonian and its eigenvalues are not printed; the coherence is what they
    # produce, so agreeing there is agreeing on both.
    for suffix in ("_noDiv", "_wiDiv"):
        pa, pb = args.a_out + suffix, args.b_out + suffix
        try:
            la = open(pa).read().split("\n")
            lb = open(pb).read().split("\n")
        except OSError as exc:
            check(False, "could not read the coherence output (%s)" % exc)
            continue
        va = [float(x) for line in la for x in FLOAT.findall(line)]
        vb = [float(x) for line in lb for x in FLOAT.findall(line)]
        if len(va) != len(vb) or not va:
            check(False, "%s : coherence files have different shapes (%d vs %d)"
                  % (suffix, len(va), len(vb)))
            continue
        worst = max(abs(x - y) for x, y in zip(va, vb))
        identical = open(pa, "rb").read() == open(pb, "rb").read()
        check(worst < args.coh_tol,
              "%s : coherence agrees over %d values (worst = %.2e, tol %.0e)%s"
              % (suffix, len(va), worst, args.coh_tol,
                 "  [and byte-identical]" if identical else ""))

    if fails:
        print("\n  %d check(s) failed" % len(fails))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
