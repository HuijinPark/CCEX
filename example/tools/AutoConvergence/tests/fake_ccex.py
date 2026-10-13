"""A stand-in for `main.out -f input.json` with a known answer, for the runner tests (tests/test_runner.py).

Reads the child input the runner writes and writes CCEX-format curves out_state<s>_{wiDiv,noDiv} (and out_rdip<r>_... for
rdip_ladder): L(t) = exp(-(t/T2)^n) times a small state noise that depends only on (configuration, state, t), so the same
configuration and state are paired across settings, with

    T2 = T0 * f(k) * (1 + A_R exp(-rdip / L_R) + A_B exp(-rbath / L_B))

f(k) a configuration factor (spread SPREAD, larger for k >= 100). Hooks (environment): FAKE_COUNT (a file that gets one line
per call), FAKE_NAN_IDS (configurations written as nan), FAKE_BIG_IDS (configurations whose first state reaches 1e5 at the middle point), FAKE_JUMP_IDS (wiDiv of the first state 5 at one point, noDiv not), FAKE_STATE_SPREAD (T2 of each state scaled by 1 + this x a normal draw), FAKE_SLOW_IDS (configurations 3 x slower), FAKE_FLOOR (states keep a random phasor of this size once decayed), FAKE_FAIL_ONCE (a marker path: the first call writes half its files
and exits 1), FAKE_RANDOM (print the random-fallback message), FAKE_NONMONO (an rdip rung whose T2 is pushed off), FAKE_PLATEAU_BELOW (below this rbath the curve stays near 0.3).
FAKE_NAN_RUNG (that rdip rung ends in nan), FAKE_SPIKE (one point at |L| = 3), FAKE_IDENTICAL (every configuration and state the same curve), FAKE_LATE_DROP (near 0.3 until t0, then gone).
pcce.cache.file (pCCE): read when it holds this run's bath file and rbath (status=hit), else written (status=stored), as CCEX
prints it; FAKE_AVAAX_SEEN (the same for the avaax file, one line per run), FAKE_STATES_SEEN (a file that gets, per state file read, configuration, state, rbath, its sha256 or 'missing' and the step).
FAKE_KINK ("B0") from rbath B0 up, the third largest rdip of each run has T2 x 1.1 (a neighbour that fails only at a large
rbath); FAKE_SHAPE ("B0:n2") from rbath B0 up, the exponent is n2 with T2 scaled to keep the area (the decay time moves, the
area does not).
nstate 0 (no bath_state_file, the ensemble) writes out_wiDiv / out_noDiv. The input keys are CCEX's since A47 (2026-10-13):
bath_state_file, Defect.axis_file, pcce.cache.file, rbath and rdip a number or [inner, outer].
"""
from math import gamma
import hashlib
import json
import os
import re
import sys
from pathlib import Path

import numpy as np

T0, N, A_R, L_R, A_B, L_B = 1.0, 1.3, 0.5, 8.0, 0.3, 6.0


def _mix(z):
    """splitmix64 of a uint64 array."""
    z = z + np.uint64(0x9E3779B97F4A7C15)
    z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return z ^ (z >> np.uint64(31))


def uniform_at(seed, t, salt):
    """A uniform number in (0, 1) per time, the same at the same time on any grid (as CCEX gives the same value at a time
    whatever nstep, 2026-10-10)."""
    k = np.round(np.asarray(t) * 1e7).astype(np.uint64)            # 0.1 ns
    z = _mix(_mix(k) ^ _mix(np.full(k.shape, (seed * 8 + salt) % 2 ** 64, dtype=np.uint64)))
    return ((z >> np.uint64(11)).astype(np.float64) + 0.5) * 2.0 ** -53


def normal_at(seed, t, salt):
    return np.sqrt(-2 * np.log(uniform_at(seed, t, 2 * salt))) * np.cos(2 * np.pi * uniform_at(seed, t, 2 * salt + 1))


def t2(k, rbath, rdip):
    spread = 0.4 if k >= 100 else 0.08
    f = 1 + spread * np.random.default_rng(1000 + k).normal()
    bump = 0.0
    if os.environ.get("FAKE_NONMONO") and abs(rdip - float(os.environ["FAKE_NONMONO"])) < 1e-3:
        bump = 0.2
    return T0 * f * (1 + A_R * np.exp(-rdip / L_R) + A_B * np.exp(-rbath / L_B) + bump)


def write(path, t, v):
    path.write_text("".join(f"{a:15g}\t{x.real:+.10f}{x.imag:+.10f}j\t\n" for a, x in zip(t, v)))


def outer(v):
    """rbath or rdip: a number, or [inner, outer]."""
    return v[1] if isinstance(v, list) else v


def main():
    inp = json.loads(Path(sys.argv[sys.argv.index("-f") + 1]).read_text())
    inp["rbath"], inp["rdip"] = outer(inp["rbath"]), outer(inp["rdip"])
    if os.environ.get("FAKE_COUNT"):
        with open(os.environ["FAKE_COUNT"], "a") as f:
            f.write(os.getcwd() + "\n")
    if os.environ.get("FAKE_RANDOM"):
        print(os.environ["FAKE_RANDOM"] if os.environ["FAKE_RANDOM"] != "1" else "The bath state is randomly generated")
    k = int(re.search(r"(\d+)$", inp["bathfile"][0]).group(1))
    nstep, dt, nstate = int(inp["nstep"]), float(inp["deltat"]), int(inp["nstate"])
    if os.environ.get("FAKE_STATES_SEEN") and nstate:    # the state files this run reads: configuration, state, rbath, content hash
        with open(os.environ["FAKE_STATES_SEEN"], "a") as f:
            for st in range(1, nstate + 1):
                p = Path(f"{inp['bath_state_file']}{st}")
                f.write(f"{k} {st} {inp['rbath']} {hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else 'missing'} {Path.cwd().parent.name}\n")
    if os.environ.get("FAKE_AVAAX_SEEN"):                # the avaax file this run reads: configuration, rbath, content hash
        p = Path(inp.get("Defect", {}).get("axis_file", "none"))
        with open(os.environ["FAKE_AVAAX_SEEN"], "a") as f:
            f.write(f"{k} {inp['rbath']} {hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else 'missing'}\n")
    cache = inp.get("pcce", {}).get("cache", {}).get("file")
    if cache:                                            # the partition cache: rbath and the bath file stand for the partition
        c = Path(cache)
        mine = f"{inp['rbath']} {hashlib.sha256(Path(inp['bathfile'][0]).read_bytes()).hexdigest()}\n"
        if c.is_file() and c.read_text() == mine:
            print("[partition-cache] status=hit lookup_s=0 search_s=0 write_s=0 file=" + str(c))
        else:
            c.write_text(mine)
            print("[partition-cache] status=stored lookup_s=0 search_s=0 write_s=0 file=" + str(c))
    t = dt * np.arange(nstep)
    rdips = [float(inp["rdip"])] + [float(r) for r in inp.get("rdip_ladder", [])]
    rb = float(inp["rbath"])
    kink = sorted(rdips, reverse=True)[2] if os.environ.get("FAKE_KINK") and rb >= float(os.environ["FAKE_KINK"]) and len(rdips) > 2 else None
    n_exp = N
    if os.environ.get("FAKE_SHAPE") and rb >= float(os.environ["FAKE_SHAPE"].split(":")[0]):
        n_exp = float(os.environ["FAKE_SHAPE"].split(":")[1])
    nan = str(k) in os.environ.get("FAKE_NAN_IDS", "").split(",")
    big = str(k) in os.environ.get("FAKE_BIG_IDS", "").split(",")
    jump = str(k) in os.environ.get("FAKE_JUMP_IDS", "").split(",")
    files = []
    for j, r in enumerate(rdips):
        prefix = "out" if j == 0 else f"out_rdip{float(np.float32(r)):g}"
        same = bool(os.environ.get("FAKE_IDENTICAL"))      # every configuration and state the same curve
        T = t2(1 if same else k, float(inp["rbath"]), r) * (1.1 if kink is not None and r == kink else 1.0)
        T *= gamma(1 + 1 / N) / gamma(1 + 1 / n_exp)            # FAKE_SHAPE keeps the area T Gamma(1 + 1/n)
        if str(k) in os.environ.get("FAKE_SLOW_IDS", "").split(","):
            T *= 3.0                                             # a configuration that decays three times more slowly
        for s in (range(1, nstate + 1) if nstate else [0]):          # nstate 0: the ensemble, out_wiDiv
            sf = Path(f"{inp['bath_state_file']}{s}") if nstate and "bath_state_file" in inp else None
            # the noise follows the state file's content when there is one, so a state computed in a later block (the state
            # test writes states s + 1 ... as 1, 2, ...) is the same state as in a run of all of them
            seed = int(hashlib.sha256(sf.read_bytes()).hexdigest()[:12], 16) + 10 ** 6 * k if sf is not None and sf.is_file() else 10 ** 6 * k + s
            seed = 1 if same else seed
            rng = np.random.default_rng(seed)
            # an ensemble curve (nstate 0) has no state noise, as in CCEX (2026-10-10)
            noise = (0.03 if nstate else 0.0) * (normal_at(seed, t, 0) + 1j * normal_at(seed, t, 1))
            Ts = T * (1 + float(os.environ["FAKE_STATE_SPREAD"]) * rng.normal()) if os.environ.get("FAKE_STATE_SPREAD") else T
            v = np.exp(-(t / abs(Ts)) ** n_exp) * (1 + noise)
            if os.environ.get("FAKE_FLOOR") and nstate:          # a single-state curve that never goes to 0: a random phasor
                v = v + float(os.environ["FAKE_FLOOR"]) * (1 - np.exp(-(t / abs(Ts)) ** n_exp)) * np.exp(2j * np.pi * uniform_at(seed, t, 7))
            if os.environ.get("FAKE_PLATEAU_BELOW") and float(inp["rbath"]) < float(os.environ["FAKE_PLATEAU_BELOW"]):
                v = 0.3 + 0.7 * v                  # too small a bath: the coherence stays near 0.3
            if os.environ.get("FAKE_LATE_DROP"):          # stays near 0.3, then drops at t0: a plateau only in a short window
                t0 = float(os.environ["FAKE_LATE_DROP"])
                v = (0.3 + 0.7 * v) * np.exp(-(np.maximum(t - t0, 0) / 0.15) ** 2)
            if os.environ.get("FAKE_MOD"):          # an ESEEM-like modulation under the envelope: "period(ms):depth"
                per, dep = map(float, os.environ["FAKE_MOD"].split(":"))
                v = v * (1 - dep * (1 - np.cos(2 * np.pi * t / per)))
            v[0] = 1.0
            if os.environ.get("FAKE_SPIKE") and j == 0:
                v[nstep - 3] = 3.0                 # one point above 1, as a CCE division by a near-zero term gives
            if nan:
                v[nstep // 2] = np.nan
            if big and s in (0, 1):                # one state of the configuration blows up (a CCE division by a near-zero term)
                v[nstep // 2] = 1e5
            if j > 0 and os.environ.get("FAKE_NAN_RUNG") and abs(r - float(os.environ["FAKE_NAN_RUNG"])) < 1e-3:
                v[-1] = np.nan                     # one ladder rung with a non-finite last point
            g = f"_state{s}" if nstate else ""
            wi = v
            if jump and s in (0, 1):               # wiDiv jumps at one point, noDiv does not (a divergent point to correct)
                wi = v.copy(); wi[nstep // 3] = 5.0
            files += [(Path(f"{prefix}{g}_wiDiv"), wi)]
            sm = inp.get("savemode")                          # like CCEX: savemode.nodiv false writes _wiDiv alone
            if not (isinstance(sm, dict) and sm.get("nodiv", True) is False):
                files += [(Path(f"{prefix}{g}_noDiv"), v)]
    marker = os.environ.get("FAKE_FAIL_ONCE")
    if marker and not Path(marker).exists():
        Path(marker).write_text("")
        for p, v in files[:len(files) // 2]:
            write(p, t, v)
        sys.exit(1)
    for p, v in files:
        write(p, t, v)


if __name__ == "__main__":
    main()
