"""Truncation judgement: a setting (smaller rdip or rbath) against a reference ((internal notes),
registered 2026-10-07; first used as (internal notes)).

Both are compared on |mu(t)| over [0, T_d] (mu: equal-weight mean over configurations of their state means; the two sets are
paired by configuration).

  1. fit lines f = exp(-(t/T2)^n) (t2.fit); pass when max |df| <= max_tol and area (mean |df|) <= area_tol;
  (fits: autoconv.t2.fit, on the upper envelope by default; a fit that fails goes to 3)
  2. fit-quality check: s = Gaussian-smoothed (width sigma) |mu_r| - |mu_ref| must stay within 2 jackknife SE + quality_floor
     of df at every point; modulation common to both settings passes, a fit that misses a real difference does not;
  With fit_only (the default since 2026-10-07, user decisions) the fit lines decide by the relative difference of the area
  under them, |A_setting / A_reference - 1| <= rel_tol (2 %), A = T2 Gamma(1 + 1/n); when either fit fails, both curves are
  judged on the raw area int_0^T_d |mu| dt (trapezoid) with the same rel_tol, route "raw-area" (user 2026-10-12; "fit-failure"
  only when that area is not usable); the max/area of rule 1, the smoothed difference and the check are reported for
  information. Without it:
  3. when the check fails: the area of s alone (a strong modulation makes the maximum follow local spikes): > area_tol fails,
     otherwise sampling-unresolved when its jackknife SE exceeds area_tol, else pass.

Defaults: max_tol 0.01, area_tol 0.005, sigma 0.5 us, quality_floor 0.001, fit floor 0.05. Times in the CCEX unit (ms).
"""
from math import gamma

import numpy as np

from . import t2 as t2mod

DEFAULTS = dict(max_tol=0.01, area_tol=0.005, sigma=0.5e-3, quality_floor=0.001, floor=t2mod.FLOOR, points="envelope",
                n_range=(t2mod.N_MIN, t2mod.N_MAX), span=None,
                fit_only=True, rel_tol=0.02, smooth=True)


def fit_area(T2, n):
    """Area under the fitted coherence exp(-(t/T2)^n) from 0 to infinity: T2 Gamma(1 + 1/n)."""
    return T2 * gamma(1 + 1 / n)


def _smooth(times, sigma, x, rows=1024):
    """Gaussian smoothing (width sigma) of x, the kernel normalized per row; built a block of rows at a time, never the N x N
    matrix (25001 points made 5 GB, review-20261010-envelope-step #3)."""
    out = np.empty(len(times))
    for i in range(0, len(times), rows):
        k = np.exp(-0.5 * ((times[i:i + rows, None] - times[None, :]) / sigma) ** 2)
        out[i:i + rows] = (k / k.sum(1, keepdims=True)) @ x
    return out


def _fitline(times, y, floor, points, n_range=(t2mod.N_MIN, t2mod.N_MAX), span=None):
    T, n, _ = t2mod.fit(times, y, floor=floor, points=points, n_range=n_range, span=span)
    return T, n, t2mod.model(times, T, n)


def judge(times, setting, reference, decay_time, **opts):
    """setting, reference: arrays configurations x time of state means (complex), paired by row. Returns a dict with the
    route ('fit' or 'smoothed'), the verdict ('pass', 'fail', 'sampling-unresolved') and every number used."""
    o = dict(DEFAULTS, **opts)
    a_all, b_all = np.asarray(setting), np.asarray(reference)
    if a_all.shape != b_all.shape or a_all.shape[0] < 2:
        raise ValueError("setting and reference need the same shape and at least 2 paired configurations")
    m = times <= decay_time * (1 + 1e-9)
    C = a_all.shape[0]

    def smooth(idx):
        a, b = np.abs(a_all[idx].mean(0)), np.abs(b_all[idx].mean(0))
        return _smooth(times, o["sigma"], a - b), a, b

    def fitparts(idx):
        a, b = np.abs(a_all[idx].mean(0)), np.abs(b_all[idx].mean(0))
        Ta, na, fa = _fitline(times, a, o["floor"], o["points"], o["n_range"], o["span"])
        Tb, nb, fb = _fitline(times, b, o["floor"], o["points"], o["n_range"], o["span"])
        return (Ta, na), (Tb, nb), fa - fb, fa, fb

    allc = np.arange(C)
    jse = lambda arrs: np.sqrt((C - 1) / C * ((np.array(arrs) - np.mean(arrs, 0)) ** 2).sum(0))
    if o["fit_only"] and not o["smooth"]:     # the runner: the smoothed difference only informs, left out (review 10-10 #3)
        a, b = np.abs(a_all.mean(0)), np.abs(b_all.mean(0))
        s = se_s = np.full(len(times), np.nan)
        out = dict(smooth_max=None, smooth_area=None, smooth_area_se=None)
    else:
        s, a, b = smooth(allc)
        jk_s = [smooth(np.delete(allc, i))[0] for i in range(C)]
        se_s = jse(jk_s)
        out = dict(smooth_max=float(np.abs(s[m]).max()), smooth_area=float(np.abs(s[m]).mean()),
                   smooth_area_se=float(jse([np.abs(x[m]).mean() for x in jk_s])))
    out.update(decay_time=float(decay_time), points=int(m.sum()), configurations=C, settings=o)
    try:                     # a fit that fails (no decay shape to fit) sends the judgement to the smoothed route
        (Ta, na), (Tb, nb), df, fa, fb = fitparts(allc)
        jk = [fitparts(np.delete(allc, i)) for i in range(C)]
    except (ValueError, RuntimeError) as e:
        out.update(T2=None, n=None, T2_ref=None, n_ref=None, dT2_rel=None, dT2_rel_se=None, fit_max=None, fit_area=None, area_rel=None,
                   quality_ok=False, fit_error=str(e))
        df = fa = fb = np.full(len(times), np.nan)
    else:
        out.update(T2=Ta, n=na, T2_ref=Tb, n_ref=nb, dT2_rel=Ta / Tb - 1, dT2_rel_se=float(jse([x[0][0] / x[1][0] - 1 for x in jk])),
                   area_rel=fit_area(Ta, na) / fit_area(Tb, nb) - 1,
                   fit_max=float(np.abs(df[m]).max()), fit_area=float(np.abs(df[m]).mean()),
                   quality_ok=None if out["smooth_max"] is None else bool((np.abs(s - df)[m] <= 2 * se_s[m] + o["quality_floor"]).all()))
    if o["fit_only"]:        # 2026-10-07 user decisions: the T2 fit lines decide (a strong modulation spoils raw differences),
        # by the relative difference of the area under them, T2 Gamma(1 + 1/n) (T2 and n together, no blow-up near zero)
        if out["fit_max"] is None:
            # a fit failed on either side (user 2026-10-12, plan 2.1 "맞춤이 실패한 비교는 원곡선 넓이로"): both curves are judged
            # on the raw area, the trapezoid of |mean| over [0, T_d], with the same rel_tol; not judged when that is not finite
            ra, rb = float(np.trapezoid(a[m], times[m])), float(np.trapezoid(b[m], times[m]))
            rel = ra / rb - 1 if rb > 0 else float("nan")
            if np.isfinite(rel):
                out.update(route="raw-area", verdict="pass" if abs(rel) <= o["rel_tol"] else "fail",
                           raw_area=ra, raw_area_ref=rb, raw_area_rel=rel)
            else:
                out.update(route="fit", verdict="fit-failure")
        else:
            out.update(route="fit", verdict="pass" if abs(out["area_rel"]) <= o["rel_tol"] else "fail")
    elif out["quality_ok"]:
        out["route"] = "fit"
        out["verdict"] = "pass" if out["fit_max"] <= o["max_tol"] and out["fit_area"] <= o["area_tol"] else "fail"
    else:
        out["route"] = "smoothed"
        out["verdict"] = ("fail" if out["smooth_area"] > o["area_tol"] else
                          "sampling-unresolved" if out["smooth_area_se"] > o["area_tol"] else "pass")
    # curves for the report (not part of the JSON)
    out["_curves"] = dict(abs_setting=a, abs_reference=b, fit_setting=fa, fit_reference=fb, fit_diff=df, smooth_diff=s,
                          smooth_se=se_s, mask=m)
    return out
