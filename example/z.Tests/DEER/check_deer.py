#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regression checks for the DEER bath-pulse path.

Every check is either an exact identity the code must satisfy, or a comparison against the
closed form, which at cluster order 1 is not an approximation: the bath spins do not interact,
so the echo is a product over them and can be written down.

    V(t) = prod_j cos(a_j * t * min(f, 1-f)),   a_j = 2 pi * 52.04 (1-3cos^2 theta_j) / r_j^3

f is the fraction of the total time at which the bath pi lands.  At f = 0.5 the bath pulse sits
on the echo centre and the argument reduces to the familiar a_j t / 2.

Run run_test.sh first; this reads what it produced.  Exit status is 0 only if every check
passes, so it can be wired into CI.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BATH = os.path.join(HERE, "bath_DEERtest_12spin.txt")
AVAAX = os.path.join(HERE, "avaax_DEERtest")
D_MHZ = 52.04                       # (mu0/4pi) gamma_e^2 hbar / h, two electrons [MHz nm^3]
PUMPED_AXIS = 1
TOL_EXACT = 0.0                     # identities must hold bit for bit
TOL_CLOSED = 2e-4                   # closed form vs code
TOL_REFOCUS = 1e-5                  # refocusing is exact in theory; this is the
                                    # propagator's own precision with a 2.88 GHz ZFS


def read_curve(case):
    p = os.path.join(HERE, case, "res_CCE", "test_wiDiv")
    t, v = [], []
    for line in open(p):
        w = line.split()
        if len(w) == 2:
            t.append(float(w[0]) * 1e3)                       # ms -> us
            v.append(complex(w[1].replace("+-", "-")).real)
    return np.array(t), np.array(v)


def closed_form(t, frac):
    rows = [l.split() for l in open(BATH).read().splitlines()[1:]]
    rows = [r for r in rows if len(r) == 4]
    xyz = np.array([[float(x) for x in r[:3]] for r in rows])
    ax = np.array([int(round(float(x))) for x in open(AVAAX).read().split()[1:]])
    r = np.linalg.norm(xyz, axis=1)
    sel = ax == PUMPED_AXIS
    f = D_MHZ * (1 - 3 * (xyz[sel, 2] / r[sel]) ** 2) / (r[sel] / 10.0) ** 3   # MHz
    eff = min(frac, 1.0 - frac)
    return np.prod(np.cos(2 * np.pi * np.outer(t, f) * eff), axis=1)


def report(name, ok, detail):
    print("  [{}] {:<46} {}".format("PASS" if ok else "FAIL", name, detail))
    return ok


def main():
    for c in ("nopulse", "onres", "offres", "multi_hit", "multi_miss", "nearpulse"):
        if not os.path.exists(os.path.join(HERE, c, "res_CCE", "test_wiDiv")):
            sys.exit("missing output for '{}' - run run_test.sh first".format(c))

    t, none = read_curve("nopulse")
    _, on = read_curve("onres")
    _, off = read_curve("offres")
    _, mhit = read_curve("multi_hit")
    _, mmiss = read_curve("multi_miss")
    _, near = read_curve("nearpulse")
    ok = []

    print("\nDEER bath-pulse regression, 12 spins, gCCE order 1, -N 0\n")

    ok.append(report("order 1 without a pump refocuses exactly",
                     np.abs(none - 1.0).max() <= TOL_REFOCUS,
                     "max |V - 1| = {:.2e}".format(np.abs(none - 1.0).max())))

    ok.append(report("off-resonant pump is identical to no pump",
                     np.abs(off - none).max() <= TOL_EXACT,
                     "max |diff| = {:.2e}".format(np.abs(off - none).max())))

    ok.append(report("multi-tone [resonant, far] == single resonant",
                     np.abs(mhit - on).max() <= TOL_EXACT,
                     "max |diff| = {:.2e}".format(np.abs(mhit - on).max())))

    ok.append(report("multi-tone [far, far] == no pump",
                     np.abs(mmiss - none).max() <= TOL_EXACT,
                     "max |diff| = {:.2e}".format(np.abs(mmiss - none).max())))

    cf = closed_form(t, 0.5)
    e = np.abs(on - cf).max()
    ok.append(report("pump at the echo centre matches the closed form",
                     e <= TOL_CLOSED, "max |diff| = {:.2e}".format(e)))

    # Regression for the bath-pulse timing tolerance.  A bath pi 0.005 of the sequence away from
    # the qubit pi used to match at both segment edges and be applied twice; two pi pulses cancel,
    # so the pump silently did nothing and this curve collapsed onto 'nopulse'.
    cfn = closed_form(t, 0.505)
    e = np.abs(near - cfn).max()
    ok.append(report("pump next to the qubit pulse is applied once, not twice",
                     e <= TOL_CLOSED, "max |diff| vs closed form = {:.2e}".format(e)))
    d = np.abs(near - none).max()
    ok.append(report("  ... and is therefore not equal to no pump",
                     d > 1e-3, "max |diff| vs nopulse = {:.2e}".format(d)))

    bad = []
    for c in ("nopulse", "onres", "offres", "multi_hit", "multi_miss", "nearpulse"):
        a = os.path.join(HERE, c, "res_CCE", "test_wiDiv")
        b = a.replace("_wiDiv", "_noDiv")
        if os.path.exists(b) and open(a).read() != open(b).read():
            bad.append(c)
    ok.append(report("order 1 has no inclusion-exclusion division (wiDiv == noDiv)",
                     not bad, "differing: {}".format(bad if bad else "none")))

    n = sum(ok)
    print("\n{} / {} checks passed\n".format(n, len(ok)))
    return 0 if n == len(ok) else 1


if __name__ == "__main__":
    sys.exit(main())
