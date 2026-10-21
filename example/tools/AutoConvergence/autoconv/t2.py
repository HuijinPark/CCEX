"""T2: the reported quantity the sampling step is judged on (2026-10-07 user decision, (internal notes)).

T2 and n come from fitting exp(-(t/T2)^n) to |mu(t)|, mu the pooled complex mean, on the points from t = 0 until |mu| first falls
below `floor` (0.05). The sampling error of T2 is estimated by linearization: each state curve L_cs moves T2 by
psi_cs = dT2 / d(mu) . (L_cs - mu) (central difference, the fitted points fixed), and the one-way layout of psi gives

    W_T = mean_c var_s psi_cs          within-configuration
    B_T = max(var_c mean_s psi_cs - W_T / S, 0)
    SE(C, S) / T2 = sqrt(B_T / C + W_T / (C S)) / T2

A two-stage bootstrap (configurations, then states within them) of the present sample checks it; on the CR0 pilot the
linearization was about 10 % below the bootstrap (8.6 % vs 9.5 %), which is why the bootstrap is reported next to it.
Times are in the unit of the CCEX output (ms); T2 is returned in the same unit.

Since 2026-10-07 (user decision) the quantity is the area under the fit, A = T2 Gamma(1 + 1/n) (the one the truncation steps
compare): psi, B, W, the SE and the bootstrap are of A, and the SE is relative to A; T2 and n are kept for the report.

The fit points (user 2026-10-10, plan 2.1): points="peaks" (the runner's) fits t = 0 and the points higher than both
neighbours, chosen on each curve; fewer than 3 peaks fit every point. The T2 is a measure for the convergence checks, not the
exact envelope. points="envelope" (the analyze driver) steps through upper_points.
"""
import math
import warnings

import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit

FLOOR = 0.05
N_MIN, N_MAX = 0.3, 4.0


def model(t, T, n):
    return np.exp(-(t / T) ** n)


def fit_points(y, floor=FLOOR):
    """Number of points from t = 0 until |y| first falls below floor (all points if it never does)."""
    below = np.nonzero(np.abs(y) < floor)[0]
    return int(below[0]) if len(below) else len(y)


def upper_points(times, y, space=3):
    """Indices of the upper envelope of a modulated curve: from each point, step to the point among the next `space` with
    the largest slope (the least decrease). Ported from cce_analysis.t2._pick_max_slope (~/bin/code_develop/cce-analysis,
    itself cce_fit.pick_high_data); used for curves with a strong modulation (NV-13C ESEEM, 2026-10-07 user decision)."""
    n = len(times)
    idx, i = [0], 0
    while i < n - 1:
        best_j, best = i + 1, (y[i + 1] - y[i]) / (times[i + 1] - times[i])
        for j in range(i + 1, min(i + 1 + space, n)):
            d = (y[j] - y[i]) / (times[j] - times[i])
            if d > best:
                best, best_j = d, j
        i = best_j
        idx.append(i)
    return np.array(sorted(set(idx)))


OUT_REL = 1e-5     # relative rounding of the times CCEX writes (6 significant digits), with a margin


def steps(times, span):
    """Points ahead for a time span h (ms) on this grid: floor(h / dt), dt from the first and last time. CCEX writes times
    with 6 significant digits (%15g), so a ratio within OUT_REL of a whole number counts as that number (review
    20261010-160920 R1: the same deltat on 201 points read as a step longer than the span from 401 points). Fewer than 1:
    the grid is coarser than the span, a ValueError."""
    dt = (float(times[-1]) - float(times[0])) / (len(times) - 1)
    r = span / dt
    s = int(math.floor(r * (1 + OUT_REL)))
    if s < 1:
        raise ValueError(f"the time step {dt:.4g} ms is longer than the envelope span {span:.4g} ms")
    return s


def envelope_idx(times, y, floor, space):
    """The envelope points up to the first one below floor (all of them if none is)."""
    up = upper_points(times, y, space)
    below = np.nonzero(y[up] < floor)[0]
    return up[:below[0]] if len(below) else up


MIN_PEAKS = 3              # fewer peaks than this before the floor: the curve has no modulation to follow, every point is fitted


def peaks(y, floor=FLOOR):
    """t = 0 and the points higher than both neighbours (y[i] > y[i-1], y[i] >= y[i+1]) up to the first one below floor
    (user 2026-10-10)."""
    i = np.arange(1, len(y) - 1)
    p = np.r_[0, i[(y[i] > y[i - 1]) & (y[i] >= y[i + 1])]]
    below = np.nonzero(y[p] < floor)[0]
    p = p[:below[0]] if len(below) else p
    return p


def peak_idx(y, floor=FLOOR):
    """The fit points: peaks(), or every point up to the floor when fewer than MIN_PEAKS peaks."""
    p = peaks(y, floor)
    return p if len(p) - 1 >= MIN_PEAKS else np.arange(fit_points(y, floor))


def fit(times, mu, npts=None, floor=FLOOR, points="envelope", idx=None, n_range=(N_MIN, N_MAX), span=None):
    """(T2, n, idx) of exp(-(t/T2)^n) fitted to |mu|. points="envelope": the upper envelope (upper_points) up to the first
    envelope point below floor; "all": every point up to the first point below floor. idx (from an earlier call) fixes the
    fitted points, e.g. when the curve is perturbed for the linearization; npts (an int) is the old form of idx = 0..npts-1.
    span: the envelope step as a time (ms, steps()); None keeps the old 3 points (the analyze driver)."""
    y = np.abs(mu)
    if idx is None and npts is not None:
        idx = np.arange(npts)
    if idx is None:
        if points == "all":
            idx = np.arange(fit_points(y, floor))
        elif points == "peaks":
            idx = peak_idx(y, floor)
        elif points == "envelope":
            idx = envelope_idx(times, y, floor, 3 if span is None else steps(times, span))
        else:
            raise ValueError("points must be 'peaks', 'envelope' or 'all'")
    if len(idx) < 3:
        raise ValueError("fewer than 3 points before |mu| falls below the fit floor")
    e = np.nonzero(y[idx] < math.exp(-1))[0]
    t0 = times[idx][e[0]] if len(e) else times[idx][-1]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", OptimizeWarning)     # an exact fit has no covariance; only p is used
        warnings.simplefilter("ignore", RuntimeWarning)      # trial parameters of the optimizer may overflow
        p, _ = curve_fit(model, times[idx], y[idx], p0=(t0, 1.2), maxfev=20000)
    # T2 > 0 and a stretching exponent in [N_MIN, N_MAX], finite (review 20261007i: the fit had no bounds); a fit outside
    # is a failure, not a value (bounding the optimizer instead would change the converged answers recorded so far).
    lo, hi = n_range                 # t2_fit.n_range of the auto_convergence block (2026-10-08)
    if not (np.isfinite(p).all() and p[0] > 0 and lo <= p[1] <= hi):
        raise ValueError(f"T2 fit outside its range (T2 {p[0]:.4g}, n {p[1]:.4g}; T2 > 0, {lo:g} <= n <= {hi:g})")
    return float(p[0]), float(p[1]), idx


def components(times, data, eps=1e-2, floor=FLOOR, points="envelope", quantity="area", n_range=(N_MIN, N_MAX), span=None):
    """Linearized variance components of the area A (quantity "T2": of T2, as recorded before 2026-10-07) from the balanced
    sample (configurations x states x time). se_rel is relative to comp["A"], which holds T2 itself for quantity "T2"."""
    x = np.array([d for d in data])
    c, s = x.shape[:2]
    if c < 2:
        raise ValueError("T2 sampling needs at least 2 configurations")
    mu = x.mean((0, 1))
    T2, n, idx = fit(times, mu, floor=floor, points=points, n_range=n_range, span=span)
    area = (lambda y: _area(*fit(times, y, idx=idx, n_range=n_range)[:2])) if quantity == "area" else (lambda y: fit(times, y, idx=idx, n_range=n_range)[0])
    psi = np.empty((c, s))
    for i in range(c):
        for j in range(s):
            d = x[i, j] - mu
            psi[i, j] = (area(mu + eps * d) - area(mu - eps * d)) / (2 * eps)
    w = float(np.mean(psi.var(1, ddof=1))) if s > 1 else 0.0     # one curve per configuration (ensemble): B holds it all
    b_raw = float(psi.mean(1).var(ddof=1)) - w / s
    return dict(T2=T2, n=n, A=_area(T2, n) if quantity == "area" else T2, quantity=quantity, fit_points=len(idx), points=points, B=max(b_raw, 0.0), W=w, negative_B=b_raw < 0,
                C=c, S=s)


def _area(T2, n):
    """Area under exp(-(t/T2)^n) from 0 to infinity (truncation.fit_area)."""
    return T2 * math.gamma(1 + 1 / n)


def se_rel(comp, c, s):
    """Relative SE of the area A with c configurations of s states."""
    return math.sqrt(comp["B"] / c + comp["W"] / (c * s)) / comp["A"]


MAX_BOOT_FAILED = 0.05      # bootstrap resamples whose fit fails: left out up to this share, more is a failure (user 2026-10-10)


def _boot_areas(times, x, draw, n_boot, rng, floor, points, n_range, info, span=None):
    """Areas of n_boot resamples (draw(rng) -> the resampled curves); a resample whose fit fails is left out and counted
    (P1 CR0 cCCE2 single, idm 17061: one of 400 left 2 fit points and stopped the step); more than MAX_BOOT_FAILED raises."""
    boot, failed = [], 0
    for _ in range(n_boot):
        y = draw(rng)
        try:
            boot.append(_area(*fit(times, y.mean((0, 1)), floor=floor, points=points, n_range=n_range, span=span)[:2]))
        except (ValueError, RuntimeError):
            failed += 1
    if info is not None:
        info.update(bootstrap_failed=failed, bootstrap_n=n_boot)
    if failed > MAX_BOOT_FAILED * n_boot or len(boot) < 2:
        raise ValueError(f"{failed} of {n_boot} bootstrap resamples could not be fitted (more than {MAX_BOOT_FAILED:.0%})")
    return boot


def bootstrap_se_rel(times, data, n_boot=400, seed=20261007, floor=FLOOR, points="envelope", n_range=(N_MIN, N_MAX), info=None, span=None):
    """Relative SE of the area A from a two-stage bootstrap of the present sample (configurations, then states). info (a
    dict) gets bootstrap_failed, the resamples left out because their fit failed."""
    x = np.array([d for d in data])
    c, s = x.shape[:2]
    rng = np.random.default_rng(seed)
    A = _area(*fit(times, x.mean((0, 1)), floor=floor, points=points, n_range=n_range, span=span)[:2])
    draw = lambda r: np.array([x[k, r.integers(0, s, s)] for k in r.integers(0, c, c)])
    boot = _boot_areas(times, x, draw, n_boot, rng, floor, points, n_range, info, span)
    return float(np.std(boot, ddof=1)) / A


def se_rel_states(comp):
    """The states' share of the relative SE of A: the configurations kept, only the states drawn again (2026-10-10)."""
    return math.sqrt(comp["W"] / (comp["C"] * comp["S"])) / comp["A"]


def bootstrap_se_rel_states(times, data, n_boot=400, seed=20261007, floor=FLOOR, points="envelope", n_range=(N_MIN, N_MAX),
                            info=None, span=None):
    """The states' share of the relative SE of A by the bootstrap: every configuration kept, its states drawn again."""
    x = np.array([d for d in data])
    c, s = x.shape[:2]
    rng = np.random.default_rng(seed)
    A = _area(*fit(times, x.mean((0, 1)), floor=floor, points=points, n_range=n_range, span=span)[:2])
    draw = lambda r: np.array([x[k, r.integers(0, s, s)] for k in range(c)])
    boot = _boot_areas(times, x, draw, n_boot, rng, floor, points, n_range, info, span)
    return float(np.std(boot, ddof=1)) / A


def allocate(comp, target, cost_config, cost_state, s_candidates=(2, 3, 5, 10, 20, 30, 60), max_configurations=None,
             max_states=None, correction=1.0, min_configurations=2):
    """Cheapest (C, S) with correction * SE(C, S)/T2 <= target and C <= max_configurations; cost C (a + b S). correction
    (>= 1) scales the linearized SE, e.g. by bootstrap / linearized at the pilot."""
    if not (math.isfinite(cost_config) and math.isfinite(cost_state)) or cost_config < 0 or cost_state <= 0:
        raise ValueError("costs must be finite, non-negative (configuration) and positive (state)")
    if max_states is not None:
        s_candidates = tuple(s for s in s_candidates if s <= max_states)
    table = []
    for s in s_candidates:
        need = (comp["B"] + comp["W"] / s) * (correction / (target * comp["A"])) ** 2
        c = max(min_configurations, math.ceil(need - 1e-9))   # never fewer than min_configurations (user 2026-10-07)
        table.append(dict(S=s, C=c, cost=c * (cost_config + cost_state * s), se_rel=correction * se_rel(comp, c, s),
                          feasible=max_configurations is None or c <= max_configurations))
    ok = [r for r in table if r["feasible"]]
    if ok:
        return dict(status="feasible", best=min(ok, key=lambda r: (r["cost"], r["C"] * r["S"])), table=table,
                    correction=correction, max_configurations=max_configurations)
    if not table or max_configurations is None or max_configurations < 2:
        return dict(status="source-limited", best=None, table=table, correction=correction, max_configurations=max_configurations)
    c = max_configurations
    best = min((dict(S=s, C=c, cost=c * (cost_config + cost_state * s), se_rel=correction * se_rel(comp, c, s)) for s in s_candidates),
               key=lambda r: (r["se_rel"], r["cost"]))
    return dict(status="infeasible: the source has too few configurations for the target", best=best, table=table,
                correction=correction, max_configurations=max_configurations)
