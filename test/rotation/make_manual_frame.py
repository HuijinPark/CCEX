#!/usr/bin/env python3
"""
Build the "A" side of the A/B equivalence test: the same physical system, written
out by hand in the qubit-aligned frame, so that CCEX can run it with no rotation
option at all.

    A = hand-transformed input, rotation off
    B = original input        , rotation on

If the two agree, then nothing frame-dependent was missed -- which is the one thing
the component-by-component checks cannot tell you.

What gets transformed, and about which point:

  bath file        positions are ABSOLUTE, so r -> r0 + R(r - r0)
  qubit file       r0 is the fixed point, so it is copied unchanged
  HF / QD file     the position columns are RELATIVE to the qubit
                   (reader.cpp builds rxyz = spxyz - qxyz before the lookup),
                   so they rotate about the origin: v -> R v
                   the 3x3 tensor blocks -> R T R^T
                   the "vN :" boundary vertices are relative too -> R v
                   HF "iso" is a scalar (Fermi contact) and does not move;
                   the QD "etc____" fallback rows carry a full tensor and do

Everything else -- headers, g-factor / eQ rows, atom names -- is copied verbatim.
Numbers are written with 17 significant digits so that the double CCEX parses back
is the one computed here, not a rounded version of it.

  ./make_manual_frame.py --bath-axis 0 0 1 --qubit-axis 1 1 1 --r0 X Y Z \\
        --bath in out --tensor 13 bath in out --tensor 12 bath in out
"""

import argparse
import math
import re
import sys

VERTEX = re.compile(r"^(v[1-8]\s*:\s*)(.*)$")
NUM = r"-?\d+\.?\d*(?:[eE][-+]?\d+)?"


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


def rot_vec(R, v):
    return [sum(R[i][j] * v[j] for j in range(3)) for i in range(3)]


def rot_tensor(R, T9):
    M = [[T9[3 * i + j] for j in range(3)] for i in range(3)]
    return [sum(R[i][k] * M[k][l] * R[j][l] for k in range(3) for l in range(3))
            for i in range(3) for j in range(3)]


def fmt(x):
    return "%.17g" % x


def row(name, values):
    """One tensor-file data row.

    Space separated, deliberately. READ_Tensor tokenizes with strtok(line, " ") --
    spaces only, never tabs. The original files get away with tab separators because
    every number is also space-padded; a purely tab-separated line comes back as one
    giant token, the column count never matches, and the row is dropped without a
    word. The failure then surfaces much later, as "the nuclear spin is within the
    range, but doesn't exist in A-file".
    """
    return name + "".join("  " + fmt(x) for x in values)


def transform_bath(R, r0, src, dst):
    """First line is the count header; the rest are 'x y z name'."""
    out = []
    with open(src) as fh:
        lines = fh.read().splitlines()
    out.append(lines[0])
    n = 0
    for line in lines[1:]:
        parts = line.split()
        if len(parts) != 4:
            out.append(line)
            continue
        v = [float(parts[i]) - r0[i] for i in range(3)]
        w = rot_vec(R, v)
        out.append("%s\t%s\t%s\t%s" % (fmt(w[0] + r0[0]), fmt(w[1] + r0[1]),
                                       fmt(w[2] + r0[2]), parts[3]))
        n += 1
    open(dst, "w").write("\n".join(out) + "\n")
    return n


def transform_tensorfile(R, ncol, src, dst, frame="bath"):
    """ncol = 13 for hyperfine (x,y,z,iso,+9), 12 for quadrupole (x,y,z,+9).

    frame says what the A side needs the COMPONENTS to be:
      "bath"   the file is declared bath-frame, so CCEX would transform it -- do the
               same here, R T R^T
      "qubit"  the file is declared already-computational, so CCEX leaves it alone --
               leave it alone here too

    The POSITION columns and the boundary vertices rotate either way: they are what
    the lookup matches against, and on the A side the bath has already moved.
    """
    rot_T = (lambda T: rot_tensor(R, T)) if frame == "bath" else (lambda T: T)
    out, nrow, nvert, netc = [], 0, 0, 0
    with open(src) as fh:
        lines = fh.read().splitlines()

    for line in lines:
        m = VERTEX.match(line)
        if m:
            v = [float(x) for x in m.group(2).split()]
            out.append(m.group(1) + " ".join(fmt(c) for c in rot_vec(R, v)))
            nvert += 1
            continue

        parts = line.split("\t")
        nums = [p for p in parts if re.fullmatch(r"\s*" + NUM + r"\s*", p)]

        # A data row: a name followed by exactly ncol numbers.
        if len(nums) == ncol and "etc____" not in line and "g-factor" not in line and "eQ___" not in line:
            vals = [float(p) for p in nums]
            xyz = rot_vec(R, vals[0:3])
            if ncol == 13:
                tail = [vals[3]] + rot_T(vals[4:13])   # iso is a scalar, no frame
            else:
                tail = rot_T(vals[3:12])
            out.append(row(parts[0], xyz + tail))
            nrow += 1
            continue

        # A quadrupole "etc____" fallback row carries a full tensor; the hyperfine
        # one carries a single isotropic number and does not move.
        if "etc____" in line and len(nums) == 9:
            rotated = rot_T([float(p) for p in nums])
            head = line.split("---")[0] + "---       ---       ---     "
            out.append(row(head, rotated))
            netc += 1
            continue

        out.append(line)

    open(dst, "w").write("\n".join(out) + "\n")
    return nrow, nvert, netc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bath-axis", nargs=3, type=float, required=True)
    ap.add_argument("--qubit-axis", nargs=3, type=float, required=True)
    ap.add_argument("--r0", nargs=3, type=float, required=True)
    ap.add_argument("--bath", nargs=2, action="append", default=[],
                    metavar=("SRC", "DST"))
    ap.add_argument("--tensor", nargs=4, action="append", default=[],
                    metavar=("NCOL", "FRAME", "SRC", "DST"))
    args = ap.parse_args()

    R = build_R(args.bath_axis, args.qubit_axis)
    print("  R rows: %s" % [[round(c, 9) for c in row] for row in R])

    for src, dst in args.bath:
        n = transform_bath(R, args.r0, src, dst)
        print("  bath   %-28s -> %-28s  %d spins" % (src.split("/")[-1], dst.split("/")[-1], n))

    for ncol, frame, src, dst in args.tensor:
        if frame not in ("bath", "qubit"):
            print("  ERROR: --tensor FRAME must be bath or qubit, got %r" % frame)
            return 1
        nrow, nvert, netc = transform_tensorfile(R, int(ncol), src, dst, frame)
        print("  tensor %-24s -> %-24s  %d rows, %d vertices, %d etc  (components: %s)"
              % (src.split("/")[-1], dst.split("/")[-1], nrow, nvert, netc,
                 "rotated" if frame == "bath" else "kept"))
        if nrow == 0:
            print("  ERROR: no %s-column data rows matched in %s" % (ncol, src))
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
