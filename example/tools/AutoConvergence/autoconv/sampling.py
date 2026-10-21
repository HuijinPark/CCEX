"""Sampling sizes: how many bath configurations C and states per configuration S (plan section 0.3 and 0.6).

The target is the pooled complex mean mu(t) = (1/C) sum_c (1/S_c) sum_s L_cs(t) (configurations weighted equally). Per stored
time point, method of moments for a balanced one-way layout with complex values:

    W(t) = sum_cs |L_cs - m_c|^2 / (C (S - 1))          within-configuration (state) variance
    V(t) = sum_c  |m_c  - m|^2   / (C - 1)              variance of the configuration means
    B(t) = max(V - W/S, 0)                              between-configuration variance (negative estimates cut at 0)

The curve-level error is summarised by

    sampling_rms_proxy(C, S) = sum_j w_j sqrt(B_j/C + W_j/(C S))      (area: w_j = 1/T; max: the largest term)

an expected-error scale (by Jensen, an upper bound on the expected area error when the variances are the true ones). It is
NOT a standard error or a confidence limit of the area norm; `bootstrap_norm_error` gives the resampled distribution of the
norm error of the current mean, which is what to quote for the present sample. Unbalanced inputs are reduced to their
common number of states for the moments.
"""
import math

import numpy as np


def _balanced(data):
    """config x state x time array with the common number of states (the first S of every configuration)."""
    s = min(x.shape[0] for x in data)
    if s < 2:
        raise ValueError("every configuration needs at least 2 states to estimate the within-configuration variance")
    if len(data) < 2:
        raise ValueError("at least 2 configurations are needed to estimate the between-configuration variance")
    return np.array([x[:s] for x in data])


def variance_components(data):
    """Per time point W, B (cut at 0), V, the number of negative B estimates, C and S."""
    x = _balanced(data)
    c, s, _ = x.shape
    mc = x.mean(1)
    m = mc.mean(0)
    w = (np.abs(x - mc[:, None]) ** 2).sum((0, 1)) / (c * (s - 1))
    v = (np.abs(mc - m) ** 2).sum(0) / (c - 1)
    b_raw = v - w / s
    return dict(W=w, B=np.maximum(b_raw, 0.0), V=v, negative_B=int((b_raw < 0).sum()), C=c, S=s)


def pooled_within(data):
    """Within-configuration variance per time point from every state of every configuration (unequal state counts
    allowed): sum_c sum_s |L_cs - m_c|^2 / sum_c (S_c - 1). Used to check the balanced estimate with the extra states a
    pilot adds to a few configurations (they do not enter the balanced components)."""
    num = sum((np.abs(x - x.mean(0)) ** 2).sum(0) for x in data if x.shape[0] >= 2)
    dof = sum(x.shape[0] - 1 for x in data if x.shape[0] >= 2)
    if dof == 0:
        raise ValueError("no configuration has 2 or more states")
    return num / dof


def _reduce(values, norm, mask=None):
    """area: mean over the evaluated time points; max: their maximum. mask selects the evaluated points (the decay part
    [0, T_d]; 2026-10-07 decision) so that a long near-zero tail does not dilute the area error; None = all points."""
    v = values if mask is None else values[mask]
    return float(v.mean()) if norm == "area" else float(v.max())


def sampling_rms_proxy(comp, c, s, norm="area", mask=None):
    return _reduce(np.sqrt(comp["B"] / c + comp["W"] / (c * s)), norm, mask)


def configurations_needed(comp, s, target, norm="area", mask=None):
    """Smallest C with sampling_rms_proxy(C, s) <= target. For area the proxy is (1/sqrt C) sum_j w_j sqrt(B_j + W_j/s);
    for max it is (1/sqrt C) max_j sqrt(...); both scale as 1/sqrt(C)."""
    unit = _reduce(np.sqrt(comp["B"] + comp["W"] / s), norm, mask)
    need = (unit / target) ** 2
    if not math.isfinite(need) or need > 1e15:
        raise OverflowError("the target needs an unrepresentable number of configurations")
    return max(2, math.ceil(need))


def allocate(comp, target, cost_config, cost_state, s_candidates=(2, 3, 5, 10, 20, 30, 60), norm="area",
             max_configurations=None, max_states=None, mask=None):
    """Cheapest (C, S) among the candidate S that meets the target, with cost C (a + b S); a = cost of preparing a
    configuration (generation, reading, clustering), b = cost of one state. S starts at 2 because the analysis needs two
    states per configuration to separate the variance components. max_configurations is the number of independent
    configurations the source can give (a fixed file list; None for a generator): allocations above it are infeasible,
    and when no candidate is feasible the best feasible one is reported with its predicted shortfall instead.
    The continuous optimum for the area surrogate sum(B)/C + sum(W)/(C S) is S* = sqrt(a W / (b B)) (reported for reference)."""
    if not (math.isfinite(cost_config) and math.isfinite(cost_state)) or cost_config < 0 or cost_state <= 0:
        raise ValueError("costs must be finite, non-negative (configuration) and positive (state)")
    if max_states is not None:            # the state source has only so many reproducible states (review 20261007f F4)
        s_candidates = tuple(s for s in s_candidates if s <= max_states)
        if not s_candidates:
            return dict(status="source-limited: fewer than 2 reproducible states per configuration", best=None, table=[],
                        max_configurations=max_configurations, max_states=max_states, S_star_continuous=None, S_star_note="")
    if max_configurations is not None and max_configurations < 2:
        return dict(status="source-limited: fewer than 2 configurations available, no allocation can be evaluated",
                    best=None, table=[], max_configurations=max_configurations, S_star_continuous=None, S_star_note="")
    table = []
    for s in s_candidates:
        try:
            c = configurations_needed(comp, s, target, norm, mask)
        except OverflowError:
            table.append(dict(S=s, C=None, cost=None, proxy=None, feasible=False))
            continue
        feasible = max_configurations is None or c <= max_configurations
        table.append(dict(S=s, C=c, cost=c * (cost_config + cost_state * s), proxy=sampling_rms_proxy(comp, c, s, norm, mask),
                          feasible=feasible))
    ok = [r for r in table if r["feasible"]]
    if not ok and max_configurations is None:
        return dict(status="infeasible: the target needs an unrepresentable number of configurations", best=None, table=table,
                    max_configurations=None, S_star_continuous=None, S_star_note="")
    if ok:
        best, status = min(ok, key=lambda r: (r["cost"], r["C"] * r["S"])), "feasible"
    else:
        c = max_configurations
        cands = [dict(S=s, C=c, cost=c * (cost_config + cost_state * s), proxy=sampling_rms_proxy(comp, c, s, norm, mask)) for s in s_candidates]
        best = min(cands, key=lambda r: (r["proxy"], r["cost"]))
        best["shortfall"] = best["proxy"] - target
        status = "infeasible: the source has too few configurations for the target"
    wsum = float((comp["W"] if mask is None else comp["W"][mask]).mean())
    bsum = float((comp["B"] if mask is None else comp["B"][mask]).mean())
    if bsum > 0:
        s_star, note = max(1.0, math.sqrt(cost_config * wsum / (cost_state * bsum))), "S* = max(1, sqrt(a W / (b B)))"
    else:
        s_star, note = None, "no finite optimum: estimated between-configuration variance is 0 (only C*S matters)"
    return dict(status=status, best=best, table=table, max_configurations=max_configurations, max_states=max_states,
                S_star_continuous=s_star, S_star_note=note)


def bootstrap_norm_error(data, n_boot=2000, seed=20261007, norm="area", quantile=0.95, mask=None):
    """Distribution of the norm of (resampled pooled mean - pooled mean), resampling whole configuration-mean curves
    (keeps the time correlation; state noise enters through each configuration mean). With few configurations this is a
    rough description of the present sample, not a guarantee."""
    mc = np.array([x.mean(0) for x in data])
    m = mc.mean(0)
    rng = np.random.default_rng(seed)
    errs = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, len(mc), len(mc))
        errs[i] = _reduce(np.abs(mc[idx].mean(0) - m), norm, mask)
    return dict(mean=float(errs.mean()), rms=float(np.sqrt((errs ** 2).mean())), quantile=float(np.quantile(errs, quantile)),
                q=quantile, n_boot=n_boot, configurations=len(mc))
