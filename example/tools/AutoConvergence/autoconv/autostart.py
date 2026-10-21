"""Automatic starts of the automatic convergence ((internal notes), "자동 시작값", registered 2026-10-08): the
configurations found from the bathfile pattern and split in a seeded order, and the first rbath, time window and rdip
estimated from the bath files and the gyromagnetic ratios. Every value here is a start; the runner converges from it.
"""
import hashlib
import math
import re
from pathlib import Path

import numpy as np

N_SPINS = 300            # spins inside the first rbath (a workload seed, not a convergence criterion)
K_CLOCK = 51.4           # T_seed = K / Omega; two-point fit (NV-13C T_d 2.13 ms, P1 100 ppm 500 G T_d 12.6 us), to refit after 17018
PAIR_PHASE = 0.25        # R_seed: the bath-bath dipolar phase D_ss T_seed / R^3 reaches this
HBAR_DIP = 1.054571729   # D_ab = HBAR_DIP |g_a g_b| rad A^3 / ms for g in rad/ms/G (the dipolar tensor of hamiltonian.cpp)
GAMMA_ELECTRON = -17608.597050
FLAT_RATIO = 0.2         # shortest / longest side of the coordinate box below this: a flat bath (2D, thin slab), refused
_HERE = Path(__file__).resolve().parents     # a tool on its own may sit fewer than four folders deep (review 6548a9 R3)
ISOTOPES = (_HERE[4] if len(_HERE) > 4 else Path("/nonexistent")) / "src" / "input" / "isotope_table.cpp"   # in its checkout


def find_isotope_table(binary=None):
    """The built-in isotope table (src/input/isotope_table.cpp): next to the tool in its checkout, else in the checkout of the
    binary (bin/main.out or build/bin/main.out) - a tool folder used on its own (CCEX_AUTOCONV_DIR, a shared installation)
    is not in a checkout (hard-coding check 2026-10-09). None when neither holds it."""
    places = [_HERE[4]] if len(_HERE) > 4 else []
    if binary is not None:
        places += [Path(binary).resolve().parents[i] for i in (1, 2) if len(Path(binary).resolve().parents) > i]
    for root in places:
        f = root / "src" / "input" / "isotope_table.cpp"
        if f.is_file():
            return f
    return None


def discover(pattern, root):
    """Configuration IDs k with an existing file for the bathfile pattern (literal {k} = a non-negative integer)."""
    if pattern.count("{k}") != 1 or "{k}" in str(Path(pattern).parent):
        raise ValueError("configuration.bathfile must hold {k} once, in the file name, for the configurations to be found "
                         "(or give configuration.k_range)")
    p = Path(pattern)
    folder = p.parent if p.is_absolute() else Path(root) / p.parent
    rx = re.compile(re.escape(p.name).replace(re.escape("{k}"), "(0|[1-9][0-9]*)") + r"\Z")
    ids = sorted(int(m.group(1)) for f in (folder.iterdir() if folder.is_dir() else []) if f.is_file() and (m := rx.match(f.name)))
    return ids


def order(ids, seed):
    """The seeded order: sha256 of a fixed tag, the seed and k."""
    return sorted(ids, key=lambda k: hashlib.sha256(f"ccex-autoconv-split:{seed}:{k}".encode()).hexdigest())


def species_and_xyz(path):
    rows = [line.split() for line in Path(path).read_text().splitlines()[1:] if line.strip()]
    return [r[3] if len(r) > 3 else "" for r in rows], np.array([[float(x) for x in r[:3]] for r in rows]).reshape(-1, 3)


def gyros(gyrofile, names):
    """name -> (spin, gamma rad/ms/G, source): the gyrofile first (name spin gyro), then the built-in isotope table."""
    out = {}
    if gyrofile is not None and Path(gyrofile).is_file():
        for line in Path(gyrofile).read_text().splitlines():
            f = line.split()
            if len(f) >= 3 and f[0] in names:
                out[f[0]] = (float(f[1]), float(f[2]), "gyrofile")
    missing = [n for n in names if n not in out]
    if missing and ISOTOPES.is_file():
        for m in re.finditer(r'\{"([^"]+)",\s*([-0-9.eE+]+)f?,\s*([-0-9.eE+]+)\}', ISOTOPES.read_text()):
            if m.group(1) in missing:
                out[m.group(1)] = (float(m.group(2)), float(m.group(3)), "built-in isotope table")
    return out


def count_seeds(configs, rbathcut, r_cov, props, qubit_gamma, dm, rdipcut, n_rdip=7, n_rbath=4):
    """Starts from counted spins (plan "자동 시작값 — 2D·섞인 bath까지, 단순판" S1-S3, 2026-10-08): configs is a list of
    (species names, distances to the qubit) per test configuration; props maps a species to (spin, gamma). The rbath start
    is where the mean count first reaches N_SPINS (any geometry: bulk, layers, a qubit above a layer, clusters); the clock
    sums the species' rates from their effective densities in that sphere; rdip from the strongest like-species pair."""
    notes = []
    kept = []
    for names, d in configs:
        names, d = np.asarray(names), np.asarray(d, dtype=float)
        m = (d >= rbathcut) & (d <= r_cov)
        kept.append((names[m], np.sort(d[m]), d[m]))
    every = np.unique(np.concatenate([x[1] for x in kept])) if kept else np.array([])
    mean_count = lambda r: float(np.mean([np.searchsorted(x[1], r, side="right") for x in kept]))
    counts = np.mean([np.searchsorted(x[1], every, side="right") for x in kept], axis=0) if len(every) else np.array([])
    hit = np.nonzero(counts >= N_SPINS)[0]
    if len(hit):
        b_n = float(every[hit[0]])
        up = math.ceil(b_n * 1000) / 1000                    # rounded up: the shell that reached the count stays inside
        b_start = up if up <= r_cov else math.floor(r_cov * 1000) / 1000
    else:
        b_n = None
        b_start = math.floor(r_cov * 1000) / 1000
        notes.append(f"the bath files hold {mean_count(r_cov):.0f} spins within what they cover ({r_cov:.4g} A), fewer than {N_SPINS}: "
                     f"started at {b_start:g} A (limited by the source bath)")
    vol = 4 * math.pi / 3 * (b_start ** 3 - rbathcut ** 3)
    species = sorted({n for names, _ in configs for n in names})        # every species observed in the files
    per = {}
    omega = 0.0
    for n in species:
        S, g = props[n]
        f = S * (S + 1) / 0.75
        d_qs = HBAR_DIP * abs(qubit_gamma * g) * abs(dm) * math.sqrt(f)
        d_ss = HBAR_DIP * g * g * f
        inside = float(np.mean([np.sum((x[0] == n) & (x[2] <= b_start)) for x in kept]))
        rho = inside / vol if vol > 0 else 0.0
        first = float(min(np.asarray(d, dtype=float)[np.asarray(nm) == n].min() for nm, d in configs if np.any(np.asarray(nm) == n)))
        om = rho * (d_qs ** 2 * d_ss) ** (1 / 3)
        omega += om
        per[n] = dict(spin=S, gamma=g, d_qs=d_qs, d_ss=d_ss, count=inside, rho=rho, omega=om, first_radius=first)
        if inside == 0:
            why = "beyond what the bath files cover" if first > r_cov else ("inside rbathcut" if first < rbathcut else "farther out")
            notes.append(f"species {n} has no spin within the start {b_start:g} A (first at {first:.4g} A, {why}): left out of the "
                         "clock, kept in the calculation")
        if S != 0.5:
            notes.append(f"species {n} has spin {S:g}: its coupling is scaled by spin moments, not calibrated")
    if not (math.isfinite(omega) and omega > 0):
        raise ValueError(f"no coupling clock (qubit gamma {qubit_gamma:g}, m_alpha - m_beta {dm:g}, species {species}): "
                         "give time_window.value or a deltat")
    t_seed = K_CLOCK / omega
    r_seed = max((d["d_ss"] * t_seed / PAIR_PHASE) ** (1 / 3) for d in per.values())
    top = round(min(1.25 ** 2 * r_seed, 2 * b_start), 3)
    rungs = [round(top * 0.8 ** j, 3) for j in range(n_rdip - 1, 0, -1)]
    rungs = [r for r in rungs if np.float32(r) > np.float32(rdipcut)]
    return dict(rbath_candidates=[round(b_start / 1.25 ** j, 3) for j in range(n_rbath - 1, -1, -1)], b_n=b_n, b_start=b_start,
                species=per, omega=omega, t_seed=t_seed, r_seed=r_seed, rdip_max=top, rdip_rungs=rungs, notes=notes)


def plane_coverage(xyz, q):
    """A layer (a thin slab of spins): the thinnest axis is its normal; the coverage is sqrt(a^2 + h^2) with a the in-plane
    distance from the qubit's projection to the nearest edge and h its height above the layer (plan S4)."""
    side = xyz.max(0) - xyz.min(0)
    nrm = int(np.argmin(side))
    lat = [i for i in range(3) if i != nrm]
    a = min(min(q[i] - xyz[:, i].min(), xyz[:, i].max() - q[i]) for i in lat)
    lo, hi = float(xyz[:, nrm].min()), float(xyz[:, nrm].max())
    h = max(lo - q[nrm], q[nrm] - hi, 0.0)              # to the nearest occupied layer of a slab (0 inside it)
    return float(math.sqrt(max(a, 0.0) ** 2 + h ** 2)) if a > 0 else 0.0


