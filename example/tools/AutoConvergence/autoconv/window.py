"""Time window (plan section 0.1 and the review corrections in 0.6).

The decay is judged on the pooled complex mean over configurations and states, mu(t) (not one configuration, not mean |L|).
A first crossing of |mu| below the threshold is not enough: coherence can rise again (revivals, slow tails, sampling
noise). The envelope used here is the suffix maximum

    env(t) = max_{t' >= t, t' in window} ( |mu(t')| + k * se(t') )

with se the pointwise standard error of the pooled mean (configurations treated as the independent units). T_d is the first
stored time with env(T_d) < threshold, i.e. after T_d the curve (with its uncertainty) never comes back above the threshold
inside the computed window. Outcomes:

    resolved            T_d found and the window reaches window_margin * T_d -> recommended window = window_margin * T_d
    needs-extension     T_d found but the computed window is shorter than window_margin * T_d (re-rises after the window
                        end cannot be excluded) -> compute a longer trial window
    window-unresolved   no T_d inside the window (right-censored)

The suffix maximum only looks inside the computed window; nothing here excludes a revival after it.

Isolated spikes (2026-10-07 user decision, rule 2): a CCE expansion can put a single time point far off (a strongly coupled spin
whose own coherence passes near zero at that point; NV-13C example). Before the suffix maximum, the envelope is replaced by its
running median over `median_points` neighbouring points (9 by default; at the ends only the points that exist), so one or two
points are ignored and a revival lasting several points is kept. The times where the raw envelope is above the threshold but the
median is not are reported as `ignored_spikes`. median_points = 1 turns this off.
"""
import numpy as np


def pooled_mean(data):
    """Pooled complex mean (configurations weighted equally) and its pointwise standard error from the spread of the
    configuration means (needs at least 2 configurations; with one configuration the state spread is used instead)."""
    mc = np.array([x.mean(0) for x in data])
    mu = mc.mean(0)
    if len(mc) >= 2:
        se = np.sqrt((np.abs(mc - mu) ** 2).sum(0) / (len(mc) - 1) / len(mc))
    else:
        x = data[0]
        se = np.sqrt((np.abs(x - x.mean(0)) ** 2).sum(0) / (len(x) - 1) / len(x))
    return mu, se


def running_median(x, points):
    """Median over `points` neighbouring values (odd), using only the values that exist at the ends."""
    h = points // 2
    return x.copy() if h == 0 else np.array([np.median(x[max(0, i - h):i + h + 1]) for i in range(len(x))])


def decay_time(times, mu, se, threshold, k=2.0, median_points=9):
    """First stored time after which the running median of |mu| + k se stays below the threshold to the window end (None if
    never)."""
    upper = running_median(np.abs(mu) + k * se, median_points)
    env = np.maximum.accumulate(upper[::-1])[::-1]
    below = np.nonzero(env < threshold)[0]
    return (float(times[below[0]]), int(below[0])) if len(below) else (None, None)


def choose_window(times, data, threshold=0.01, margin=1.5, k=2.0, median_points=9):
    if median_points < 1 or median_points % 2 == 0:
        raise ValueError("median_points must be a positive odd number")
    mu, se = pooled_mean(data)
    td, _ = decay_time(times, mu, se, threshold, k, median_points)
    t_end = float(times[-1])
    first_cross = np.nonzero(np.abs(mu) < threshold)[0]
    report = dict(threshold=threshold, margin=margin, k=k, window_end=t_end,
                  first_crossing=float(times[first_cross[0]]) if len(first_cross) else None,
                  max_after_first_crossing=float(np.abs(mu[first_cross[0]:]).max()) if len(first_cross) else None,
                  median_points=median_points)
    raw = np.abs(mu) + k * se
    med = running_median(raw, median_points)
    report["ignored_spikes"] = [float(x) for x in times[(raw >= threshold) & (med < threshold)]]
    per_config = []
    for x in data:
        tc, _ = decay_time(times, x.mean(0), np.zeros(len(times)), threshold, 0.0, median_points)
        per_config.append(tc)
    report["per_configuration_decay"] = per_config
    if td is None:
        report.update(status="window-unresolved", decay_time=None, recommended_window=None)
    else:
        want = margin * td
        report.update(decay_time=td, recommended_window=want,
                      status="resolved" if t_end >= want * (1 - 1e-12) else "needs-extension")
    return report


def fit_window(times, data, threshold=0.01, margin=1.2, k=2.0, median_points=9, points="envelope", floor=0.05, look=2.5,
               n_range=(0.3, 4.0), span=None, raise_floor=True):
    """Time window from the T2 fit (2026-10-07 user decision): T_d is where the fitted line exp(-(t/T2)^n) of |mu| reaches the
    threshold, T_d = T2 (ln(1/threshold))^(1/n), and the window is margin x T_d. The fit needs the decay inside the window:
    if the upper envelope of |mu| never falls below `floor`, or the fit fails, the window needs extending. As a safeguard the
    envelope rule above is evaluated too, and the times where its judged line (median of |mu| + k SE) is above the threshold
    after T_d are reported as `late_rise` (a revival or a long tail the fitted shape does not have).

    Revivals (2026-10-07 user decision): the pilot must reach look x T_d (2.5: about 6 ms for NV-13C) so that a revival can be
    seen; until then the window needs extending. A revival is judged on the mean alone (user decision, 2026-10-07: a few slow
    configurations raise mean + 2 SE but are not a revival): if the running median of |mu| comes back above the threshold
    after T_d, the window to compute is margin x (the last such time) instead of margin x T_d. The mean + 2 SE excursions stay
    a warning (late_rise). look = 0 skips the revival check (production check).

    Noise floor (user 2026-10-10): floor and threshold are raised to 2 x the median SE of the mean over the last third of the
    window when that is higher (a sample mean of single-state curves stays near it); a plateau must be above twice that.

    Fit points (user 2026-10-10, plan 2.1): points "peaks" (the runner's) - t = 0 and the points higher than both neighbours,
    every point when fewer than 3 peaks; the decay is judged on those points. raise_floor False keeps the floor given (the
    run's, once the window has set it)."""
    from . import t2 as t2mod
    mu, se = pooled_mean(data)
    t_end = float(times[-1])
    # the noise floor (user 2026-10-10, P1 CR90 single): a sample mean cannot go below about 2 SE of itself, so the decay ends
    # there - the fit floor and the threshold are raised to 2 x the median SE over the last third when that is higher
    noise = float(np.median(2 * se[2 * (len(times) // 3):])) if len(times) >= 3 else 0.0
    threshold = max(threshold, noise)
    if raise_floor:                            # the fit floor stays the run's once set (review 20261010-first-pass-impl R4)
        floor = max(floor, noise)
    report = dict(method="fit", threshold=threshold, margin=margin, window_end=t_end, median_points=median_points, noise_level=noise,
                  fit_floor=floor)
    am = running_median(np.abs(mu), median_points)
    third = len(times) // 3
    mid, tail = float(am[third:2 * third].mean()), float(am[2 * third:].mean())
    start = float(am[:max(1, third)].max())
    # the mean stopped decaying above the floor after it had clearly decayed (review 10-07 A3: a curve still near its start
    # in a short window is not a plateau, the window needs extending): a plateau candidate (too small a bath?)
    if tail > max(floor, 2 * noise) and tail >= 0.9 * mid and tail <= 0.8 * start:     # clearly above the noise floor too
        report.update(status="plateau", decay_time=None, recommended_window=None, plateau_level=tail,
                      reason=f"|mean| stays near {tail:.3g} over the last third of the window (a plateau)")
        return report
    y = np.abs(mu)
    if points == "peaks":
        up = t2mod.peaks(y, -np.inf)           # every peak (the floor cut is the fit's): has the envelope decayed?
        up = up if len(up) - 1 >= t2mod.MIN_PEAKS else np.arange(len(times))
    elif points == "envelope":
        up = t2mod.upper_points(times, y) if span is None else t2mod.upper_points(times, y, t2mod.steps(times, span))
    else:
        up = np.arange(len(times))
    if not (np.abs(mu)[up] < floor).any():
        report.update(status="needs-extension", decay_time=None, recommended_window=None,
                      reason=f"|mean| does not fall below {floor} inside the window: nothing to fit yet")
        return report
    try:
        T2, n, idx = t2mod.fit(times, mu, floor=floor, points=points, n_range=n_range, span=span)
    except (ValueError, RuntimeError) as e:
        report.update(status="window-unresolved", decay_time=None, recommended_window=None, reason=f"T2 fit: {e}")
        return report
    td = float(T2 * np.log(1 / threshold) ** (1 / n))
    judged = running_median(np.abs(mu) + k * se, median_points)
    late = times[(times > td) & (judged >= threshold)]
    mmed = running_median(np.abs(mu), median_points)
    below = np.nonzero((times >= td) & (mmed < threshold))[0]      # a revival comes back up after the mean has gone below
    rev = times[below[0]:][mmed[below[0]:] >= threshold] if len(below) else times[:0]
    look_until = look * td
    want = margin * (float(rev[-1]) if len(rev) else td)
    need = max(want, look_until)
    used = "all" if points == "peaks" and len(t2mod.peaks(y, floor)) - 1 < t2mod.MIN_PEAKS else points
    report.update(fit_points=used, fit_npoints=int(len(idx)), n_range=tuple(n_range))
    report.update(T2=float(T2), n=float(n), decay_time=td, recommended_window=want, look_until=look_until,
                  revival=bool(len(rev)), revival_end=float(rev[-1]) if len(rev) else None,
                  status="resolved" if t_end >= need * (1 - 1e-12) else "needs-extension", needed_window=need,
                  late_rise=[float(x) for x in late], late_rise_max=float(judged[times > td].max()) if (times > td).any() else None)
    return report
