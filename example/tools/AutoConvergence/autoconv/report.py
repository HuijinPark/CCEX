"""Figures (PNG) and text summaries written next to the JSON report (--report-dir). Times are ms (CCEX); figures show µs when the time shown is below 0.1 ms (_unit).

  sampling.png / sampling.txt         |mu| with its fit and T_d; SE(T2)/T2 against the number of configurations
  truncation_<name>.png / .txt        the two |mu| with their fits; fit-line and smoothed differences with the tolerances
                                      and the jackknife band; what the fits leave out (residuals)
"""
import math
from pathlib import Path

import numpy as np

from . import t2 as t2mod

US = 1e3
TIME_LABEL = "t"     # the name of the time axis; the runner sets it from the pulse sequence (2τ for a Hahn echo, user 2026-10-07)


def _unit(t_end):
    """Display unit of the time axis from the length shown (CCEX times are ms): us below 0.1 ms (dense baths such as
    h-BN decay in microseconds), else ms. Returns (factor from ms, name)."""
    return (1e3, "µs") if t_end < 0.1 else (1.0, "ms")


L_YLIM = (-0.05, 1.08)    # |L| axes with a little room around 0 and 1 (user 2026-10-10)


def _xlim(end):
    """A time axis with a little room on both sides (user 2026-10-10)."""
    return (-0.02 * end, 1.02 * end)


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # one size for every figure (user 2026-10-10): titles, axis words and numbers alike, legends a little smaller
    plt.rcParams.update({"font.size": 18, "legend.fontsize": 16, "axes.titlesize": 18, "axes.labelsize": 18,
                         "xtick.labelsize": 18, "ytick.labelsize": 18})
    return plt


def sampling(out_dir, times, balanced, window, comp, boot_se, alloc, target, name="sampling", fixed="", figure=True):
    """<name>.png: left, the pilot curves at the chosen truncation (gray: configuration means, black: mean, red: T2 fit) up to
    the time window; right, SE(A)/A (A: area under the T2 fit) against the number of configurations - single-sample (nstate
    > 0): one line per number of states considered, the chosen one thick; ensemble (nstate 0): one line - with the target
    and the allocation (the bootstrap and the configurations available are in the text, user request 2026-10-07)."""
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    plt = _plt()
    x = np.array(balanced); mu = x.mean((0, 1))
    best = alloc.get("best") if alloc else None
    ens = x.shape[1] == 1                    # one curve per configuration: the ensemble (CCEX nstate 0)
    each = lambda c, s: f"{c} configurations" + ("" if ens else f" x {s} state" + ("s" if s != 1 else ""))
    end = window.get("recommended_window") if window else None
    f, u = _unit(end if end else times[-1])
    fig, ax = plt.subplots(1, 2, figsize=(18, 7.5))
    for m in x.mean(1):
        ax[0].plot(f * times, np.abs(m), color="gray", lw=0.8, alpha=0.5)
    ax[0].plot(f * times, np.abs(mu), "k", lw=2.0, label=f"mean ({each(x.shape[0], x.shape[1])}" + (", ensemble)" if ens else ")"))
    ax[0].plot(f * times, t2mod.model(times, comp["T2"], comp["n"]), "r-", lw=2.5,
               label=f"T2 fit: T2 {f * comp['T2']:.3g} {u}, n {comp['n']:.2f}")
    ax[0].set(xlabel=f"{TIME_LABEL} ({u})", ylabel="|L|", ylim=(0, 1.05), xlim=(0, f * (end if end else times[-1])),
              title="4. pilot" + (f" ({fixed})" if fixed else ""))
    S = best["S"] if best else x.shape[1]
    cmax = alloc.get("max_configurations") if alloc else None
    cs = np.arange(2, int(1.1 * cmax) + 1 if cmax else max(101, 2 * x.shape[0]))
    corr = alloc.get("correction", 1.0) if alloc else 1.0      # the linearized SE scaled up to the bootstrap, as the allocation
    if ens:
        ax[1].plot(cs, 100 * corr * np.array([t2mod.se_rel(comp, c, 1) for c in cs]), lw=2.5, color="C0", label="ensemble (nstate 0)")
    else:
        rows = alloc["table"] if alloc else [dict(S=S)]
        cmap = plt.get_cmap("viridis")
        for i, r in enumerate(rows):      # each number of states with the configurations it needs for the target (circle)
            sv, col = r["S"], cmap(i / max(1, len(rows) - 1) * 0.9)
            ax[1].plot(cs, 100 * corr * np.array([t2mod.se_rel(comp, c, sv) for c in cs]), color=col, lw=3.0 if sv == S else 1.3,
                       label=f"{sv} state" + ("s" if sv != 1 else "") + (" <- chosen" if sv == S else ""))
            if "C" in r and r.get("feasible", True):
                ax[1].plot(r["C"], 100 * r["se_rel"], "o", color=col, ms=11, mec="k")
    knee = None
    if not ens and alloc and len(alloc["table"]) > 1:   # one row (the states fixed, 2026-10-10): nothing to compare
        # the configurations no number of states can go below, B (correction / (target A))^2; the first number of states
        # within 1.2 x of it: more states barely reduce the configurations (marked only; the choice is the cheapest)
        floor_c = max(2, math.ceil(comp["B"] * (corr / (target * comp["A"])) ** 2 - 1e-9))
        knee = next((r for r in alloc["table"] if r.get("feasible", True) and r["C"] <= 1.2 * floor_c), None)
        if knee:
            ax[1].plot(knee["C"], 100 * knee["se_rel"], "r*", ms=24, mec="k",
                       label=f"{knee['S']} states, {knee['C']} configurations:\nmore states barely help")
    if ens and best:
        ax[1].plot(best["C"], 100 * best["se_rel"], "o", color="C0", ms=11, mec="k")
    ax[1].axhline(100 * target, color="C1", ls="--", lw=2, label=f"target {100 * target:g} %")
    ax[1].plot([], [], "o", color="white", ms=11, mec="k", label="configurations needed")
    ax[1].set(xlabel="configurations", ylabel=r"sampling error of $A = \int_0^\infty$ (T2 fit) $dt$  (%)", ylim=(0, 400 * target),
              title=f"4. {each(best['C'], best['S'])}: error {100 * best['se_rel']:.2f} %" if best
              else "4. sampling: no allocation")
    ax[0].legend(loc="upper right"); ax[1].legend(loc="upper right", ncol=1 if ens else 2, fontsize=15)
    for a in ax:
        a.grid(alpha=0.3)
    if figure:
        fig.tight_layout(); fig.savefig(out / f"{name}.png", dpi=120)
    plt.close(fig)
    lines = [f"A = T2 Gamma(1 + 1/n) = {comp['A']:.6g} ms; T2 = {comp['T2']:.6g} ms, n = {comp['n']:.4f} (fit points {comp['fit_points']})",
             f"linearized variance of A: B = {comp['B']:.4g} ms^2 (between configurations), W = {comp['W']:.4g} ms^2 (between states)",
             f"SE(A)/A at the present ({comp['C']}, {comp['S']}): linearized {t2mod.se_rel(comp, comp['C'], comp['S']):.2%}, "
             f"bootstrap {boot_se:.2%}", f"target {target:.2%}"]
    td = window.get("decay_time") if window else None
    if td is not None:
        lines.insert(0, f"window: {window['status']}, T_d {td:.4g} ms, recommended {window['recommended_window']:.4g} ms")
    if alloc:
        lines.append(f"allocation ({alloc['status']}, SE scaled by {alloc['correction']:.3f}; configurations available "
                     f"{alloc.get('max_configurations')}):")
        lines += [f"  S {r['S']:3d}: C {r['C']:4d}, SE {r['se_rel']:.2%}, cost {r['cost']:.4g}{'' if r['feasible'] else ' (too many configurations)'}"
                  for r in alloc["table"]]
        if alloc.get("best"):
            lines.append(f"  chosen: C {alloc['best']['C']} x S {alloc['best']['S']}")
        if knee:
            lines.append(f"  more states barely reduce the configurations from S {knee['S']} (C {knee['C']}; no S goes below C {floor_c})")
    (out / f"{name}.txt").write_text("\n".join(lines) + "\n")


def truncation(out_dir, name, times, j, label_setting="setting", label_reference="reference"):
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    plt = _plt()
    c = j["_curves"]; t = times * US; o = j["settings"]; td = j["decay_time"] * US
    fig, ax = plt.subplots(1, 3, figsize=(18, 4.8))
    for y, f, col, lab, T, n in ((c["abs_setting"], c["fit_setting"], "C3", label_setting, j["T2"], j["n"]),
                                 (c["abs_reference"], c["fit_reference"], "C0", label_reference, j["T2_ref"], j["n_ref"])):
        ax[0].plot(t, y, color=col, lw=0.8, alpha=0.6)
        if T is not None:
            ax[0].plot(t, f, color=col, lw=1.8, label=f"{lab}: T2 {T * US:.3f} us, n {n:.3f}")
    ax[0].set(xlabel="t (us)", ylabel="|L|", title="|mean| (thin) and fit lines (thick)", ylim=(0, 1.02))
    ax[1].fill_between(t, c["smooth_diff"] - c["smooth_se"], c["smooth_diff"] + c["smooth_se"], color="C1", alpha=0.15,
                       label="smoothed ± jackknife SE")
    ax[1].plot(t, c["smooth_diff"], "C1", lw=1.2, label=f"smoothed difference: max {j['smooth_max']:.2e}, area {j['smooth_area']:.2e}")
    if j["fit_max"] is not None:
        ax[1].plot(t, c["fit_diff"], "k", lw=1.8, label=f"fit-line difference: max {j['fit_max']:.2e}, area {j['fit_area']:.2e}")
    ax[1].axhspan(-o["max_tol"], o["max_tol"], color="green", alpha=0.08, label=f"max tolerance ±{o['max_tol']:g}")
    ax[1].axhline(0, color="k", lw=0.5)
    ax[1].set(xlabel="t (us)", ylabel=f"|L|({label_setting}) − |L|({label_reference})",
              title=f"{j['route']} route: {j['verdict']} (fit check {'ok' if j['quality_ok'] else 'failed'})")
    for y, f, col, lab in ((c["abs_setting"], c["fit_setting"], "C3", label_setting), (c["abs_reference"], c["fit_reference"], "C0", label_reference)):
        ax[2].plot(t, y - f, color=col, lw=1, label=lab)
    ax[2].axhline(0, color="k", lw=0.5)
    ax[2].set(xlabel="t (us)", ylabel="|L| − fit", title="what the fits leave out (residual)")
    for a in ax:
        a.axvline(td, color="gray", ls="--", lw=0.8); a.grid(alpha=0.3); a.legend()
    fig.tight_layout(); fig.savefig(out / f"truncation_{name}.png", dpi=130); plt.close(fig)
    fitlines = ([f"T2 {j['T2'] * US:.4f} vs {j['T2_ref'] * US:.4f} us: dT2/T2 {j['dT2_rel']:+.3%} ± {j['dT2_rel_se']:.3%}; n {j['n']:.4f} vs {j['n_ref']:.4f}",
                 f"fit-line difference: max {j['fit_max']:.3e} (tolerance {o['max_tol']:g}), area {j['fit_area']:.3e} (tolerance {o['area_tol']:g})"]
                if j["fit_max"] is not None else [f"fit failed: {j.get('fit_error')}"]
                + ([f"raw area of |L| on [0, T_d]: {j['raw_area'] * US:.4g} vs {j['raw_area_ref'] * US:.4g} us, {j['raw_area_rel']:+.3%}"]
                   if j.get("route") == "raw-area" else []))
    (out / f"truncation_{name}.txt").write_text("\n".join([
        f"{label_setting} vs {label_reference}: {j['configurations']} paired configurations, [0, T_d = {td:.3f} us] ({j['points']} points)",
        *fitlines,
        f"smoothed difference (sigma {o['sigma'] * US:g} us): max {j['smooth_max']:.3e}, area {j['smooth_area']:.3e} ± {j['smooth_area_se']:.1e}",
        f"fit check: {'ok' if j['quality_ok'] else 'failed'} -> {j['route']} route: {j['verdict']}"]) + "\n")


# ---------------------------------------------------------------- per-step figures of the --auto_conv runner (linear scales)
def window_step(out_dir, times, data, w):
    """1_window.png: the pilot's mean curve with each configuration, the T2 fit and the time window chosen."""
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    plt = _plt()
    from .window import pooled_mean
    mu, _ = pooled_mean(data)
    fig, a = plt.subplots(figsize=(9, 7.5))          # one panel the size of a panel of the two-panel figures
    win = w.get("recommended_window")
    f, u = _unit(times[-1])
    for x in data:
        a.plot(f * times, np.abs(x.mean(0)), color="gray", lw=0.5, alpha=0.6)
    a.plot(f * times, np.abs(mu), "k", lw=1.0, label="mean")
    if w.get("T2") is not None:
        a.plot(f * times, t2mod.model(times, w["T2"], w["n"]), "C3", lw=2.2, label=f"T2 fit (T2 {f * w['T2']:.3g} {u}, n {w['n']:.2f})")
    if win is not None:
        a.axvline(f * win, color="0.35", ls=":", lw=2, label=f"time window = {f * win:.3g} {u}" + (" (given)" if w.get("status") == "given" else ""))
    raised = w.get("rbath_raised") or []
    rb = f" (rbath {w['rbath']:g} Å" + (", raised from a plateau" if raised else "") + ")" if w.get("rbath") is not None else ""
    a.set(xlabel=f"{TIME_LABEL} ({u})", ylabel=f"|L({TIME_LABEL})|", ylim=L_YLIM, title=f"1. time window: {w['status']}{rb}")
    a.grid(alpha=0.3); a.legend()
    fig.tight_layout(); fig.savefig(out / "1_window.png", dpi=120); plt.close(fig)


def window_fit(out_dir, times, data, w):
    """1_window_fit.png (user 2026-10-10): the test mean, each configuration, mean +- 2 SE, the fit points (the peaks, or all
    points) with their fit, and the fit on all points for reference. w: the window report."""
    from .window import pooled_mean
    data = [x for x in data if np.isfinite(x).all()]
    if not data:
        return
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    plt = _plt()
    mu, se = pooled_mean(data); y = np.abs(mu)
    floor = w.get("fit_floor", t2mod.FLOOR)
    used = w.get("fit_points")
    f, u = _unit(times[-1])
    fig, a = plt.subplots(figsize=(13, 7.5))
    for x in data:
        a.plot(f * times, np.abs(x.mean(0)), color="0.85", lw=0.6)
    if len(data) > 1:
        a.fill_between(f * times, y - 2 * se, y + 2 * se, color="C0", alpha=0.18, label="mean ± 2 SE")
    a.plot(f * times, y, "k", lw=1.1, label=f"test mean ({len(data)} configurations, {len(times)} points)")
    idx = np.arange(0)
    if w.get("T2") is not None:
        idx = t2mod.peak_idx(y, floor) if used == "peaks" else np.arange(t2mod.fit_points(y, floor))
        a.plot(f * times[idx], y[idx], "o", color="C3", ms=5, mfc="none")
        a.plot(f * times, t2mod.model(times, w["T2"], w["n"]), "C3", lw=2.3,
               label=f"fit on {used} ({len(idx)} points): T2 {f * w['T2']:.3g} {u}, n {w['n']:.2f}")
        if used == "peaks":
            try:
                T2a, na, _ = t2mod.fit(times, mu, floor=floor, points="all", n_range=w.get("n_range", (t2mod.N_MIN, t2mod.N_MAX)))
                a.plot(f * times, t2mod.model(times, T2a, na), "C2--", lw=1.5,
                       label=f"fit on all points (reference): T2 {f * T2a:.3g} {u}, n {na:.2f}")
            except (ValueError, RuntimeError):
                pass
        head = f"fit points: {used}"
    else:
        head = "stopped: " + str(w.get("reason") or w.get("status") or "")[:80]
    a.axhline(floor, color="0.4", ls=":", lw=1.2)
    a.set(xlabel=f"{TIME_LABEL} ({u})", ylabel=f"|L({TIME_LABEL})|", ylim=L_YLIM,
          xlim=_xlim(f * (min(times[-1], 3 * times[idx[-1]]) if len(idx) else times[-1])), title=f"1. {head}")
    a.legend(); a.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(out / "1_window_fit.png", dpi=110); plt.close(fig)


def window_tried(out_dir, tried):
    """1_window.png when no time window is found (user 2026-10-10): one panel per window computed, in the order tried, the
    mean |L| of the test configurations (black) and each configuration (grey), up to the longest time computed. tried:
    [(times, [states x time per configuration], rbath)]."""
    if not tried:
        return
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    plt = _plt()
    n = len(tried); cols = min(3, n); rows = (n + cols - 1) // cols
    fig, ax = plt.subplots(rows, cols, figsize=(9 * cols, 6.5 * rows), squeeze=False)
    for a, (t, data, rb) in zip(ax.flat, tried):
        f, u = _unit(t[-1])
        for x in data:
            a.plot(f * t, np.abs(np.asarray(x).mean(0)), color="0.7", lw=0.8)
        a.plot(f * t, np.abs(np.mean([np.asarray(x).mean(0) for x in data], 0)), "k", lw=2, label=f"mean of {len(data)} configurations")
        a.set(xlabel=f"{TIME_LABEL} ({u})", ylabel=f"|L({TIME_LABEL})|", ylim=L_YLIM, xlim=_xlim(f * t[-1]),
              title=f"window {f * t[-1]:.3g} {u}, {len(t)} points, rbath {rb:g} Å")
        a.grid(alpha=0.3); a.legend()
    for a in list(ax.flat)[n:]:
        a.axis("off")
    fig.suptitle("1. time window not found: the test curves of every window tried", fontsize=20)
    fig.tight_layout(); fig.savefig(out / "1_window.png", dpi=100); plt.close(fig)


def trend_step(out_dir, number, name, times, curves, rec, chosen, window_end, rel_tol, fixed=""):
    """<number>_<name>.png for rdip (2) and rbath (3): left, |mean| of every value computed (the check and raised ones included);
    right, in %, |A(x)/A(largest) - 1| for every smaller value (A: area under the T2 fit), the tolerance and the value chosen.
    rec: the record of runner.converge (vs_reference, estimate, checks)."""
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    plt = _plt()
    values = sorted(curves)
    f, u = _unit(window_end)
    cmap = plt.get_cmap("viridis")
    fig, ax = plt.subplots(1, 2, figsize=(18, 7.5))
    for i, v in enumerate(values):
        ax[0].plot(f * times, curves[v], color=cmap(i / max(1, len(values) - 1) * 0.8), lw=3.0 if v == chosen else 1.4,
                   label=f"{name} {v:g}")
    ax[0].set(xlabel=f"{TIME_LABEL} ({u})", ylabel=f"|L({TIME_LABEL})|", ylim=L_YLIM, xlim=_xlim(f * window_end),
              title=f"{number}. {name}: curves" + (f" ({fixed})" if fixed else ""))
    pts = [e for e in rec["vs_reference"] if e["rel"] is not None]
    ax[1].plot([e["value"] for e in pts], [100 * e["rel"] for e in pts], "s:", color="C0", ms=13, lw=2,
               label=f"|A({name})/A({rec['reference']:g}) - 1|")
    raw = [e for e in pts if e.get("basis") == "raw-area"]
    if raw:                                    # a T2 fit failed: judged on the raw area (user 2026-10-12)
        ax[1].plot([e["value"] for e in raw], [100 * e["rel"] for e in raw], "o", mfc="white", mec="C3", mew=2.5, ms=17, ls="none",
                   label="raw area of |L| on [0, T_d] (a T2 fit failed)")
    ax[1].axhline(100 * rel_tol, color="C1", ls="--", lw=2, label=f"tolerance {100 * rel_tol:g} %")
    if chosen is not None:
        ax[1].axvline(chosen, color="0.35", ls=":", lw=2, label=f"chosen {name} {chosen:g}")
    e = rec.get("estimate") or {}
    c = rec["checks"][-1] if rec["checks"] else e.get("check")   # the comparison that decided
    check = (f"check |A({c['lo']:g})/A({c['hi']:g}) - 1| = {100 * c['rel']:.2f} %: {c['verdict']}"
             + (" (raw area, a T2 fit failed)" if c.get("basis") == "raw-area" else "") if c and c["rel"] is not None
             else f"check {c['lo']:g} vs {c['hi']:g}: {c['verdict']}" if c else "no check")
    # x and y from the points and the tolerance line, with room so no marker sits on the frame (user 2026-10-10)
    ys = [100 * e["rel"] for e in pts] + [100 * rel_tol]
    xs = [e["value"] for e in pts] + ([chosen] if chosen is not None else [])
    for lim, vals, pad in ((ax[1].set_ylim, ys, 0.1), (ax[1].set_xlim, xs, 0.06)):
        if not vals:                           # no comparison at all: matplotlib's own range (review R2)
            continue
        lo, hi = min(vals), max(vals)
        span = (hi - lo) or abs(hi) or 1.0
        lim(lo - pad * span, hi + pad * span)
    ax[1].set(xlabel=f"{name} (Å)", ylabel=r"relative difference of $A = \int_0^\infty$ (T2 fit) $dt$  (%)",
              title=f"{number}. {name} {chosen:g} (chosen)\n{check}" if chosen is not None else f"{number}. {name} not converged\n{check}")
    for a in ax:
        a.grid(alpha=0.3); a.legend()
    if chosen is not None:                    # the chosen value in bold, not a word (user 2026-10-10)
        for txt in ax[0].get_legend().get_texts():
            if txt.get_text() == f"{name} {chosen:g}":
                txt.set_fontweight("bold")
    fig.tight_layout(); fig.savefig(out / f"{number}_{name}.png", dpi=120); plt.close(fig)


def growth(out_dir, name, times, curves, points, threshold, window_end, what, title, chosen=None, auto_y=False):
    """<name>.png (user 2026-10-10): left, the mean |L| at every count (one colour per count, a colour bar of one box per count
    with a few numbers); right, the sampling error measured at every count, the threshold and the count chosen. curves: {count: |mean| over time}; points: [(count, error)]."""
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    plt = _plt()
    import matplotlib                    # loaded by _plt
    values = sorted(curves)
    f, u = _unit(window_end)
    n = len(values)
    cmap = matplotlib.colors.ListedColormap(plt.get_cmap("viridis")(np.linspace(0, 0.8, n)))   # no pale yellow (user 2026-10-10)
    norm = matplotlib.colors.BoundaryNorm(np.arange(n + 1) - 0.5, n)       # box i is values[i]
    fig, ax = plt.subplots(1, 2, figsize=(18, 7.5))
    for i, v in enumerate(values):
        ax[0].plot(f * times, curves[v], color=cmap(i), lw=3.0 if i == n - 1 else 1.4)
    bar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax[0], pad=0.02)
    keep = sorted({round(j * (n - 1) / 5) for j in range(6)}) if n > 6 else list(range(n))     # about 6 numbers
    bar.set_ticks(keep)
    bar.set_ticklabels([f"{values[j]:g}" for j in keep])
    bar.set_label(what)
    ax[0].set(xlabel=f"{TIME_LABEL} ({u})", ylabel=f"|L({TIME_LABEL})|", ylim=L_YLIM, xlim=_xlim(f * window_end))
    ax[0].set_title(title)
    pts = [(c, e) for c, e in points if e is not None]
    ax[1].plot([c for c, _ in pts], [100 * e for _, e in pts], "o-", color="k", ms=10, lw=2, label="measured")
    ax[1].axhline(100 * threshold, color="C1", ls="--", lw=2, label=f"threshold {100 * threshold:g} %")
    ax[1].set(xlabel=what)
    if not auto_y:                           # 4_states: the range of its points (user 2026-10-10)
        ax[1].set(xlim=(0, None), ylim=(0, None))
    ax[1].set_ylabel(r"$\sigma_A\,/\,A$  (%)")
    ax[1].set_title(f"{title.split(':')[0]}: error\n" r"$A=\int_0^\infty e^{-(t/T_2)^n}dt$ (area under the T2 fit)" "\n"
                    r"$\sigma_A$ = standard deviation of $A$")
    if chosen is not None:                   # what the step chose (user 2026-10-10)
        ax[1].axvline(chosen, color="0.35", ls=":", lw=2, label=f"chosen {chosen:g}")
    ax[1].xaxis.get_major_locator().set_params(integer=True)
    ax[0].grid(alpha=0.3); ax[1].grid(alpha=0.3); ax[1].legend(loc="best")
    fig.tight_layout(); fig.savefig(out / f"{name}.png", dpi=120); plt.close(fig)


def production_step(out_dir, times, data, T2, n, ensemble=False, fixed=""):
    """6_converged_coherence.png: the production mean curve with its T2 fit."""
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    plt = _plt()
    mu = np.mean([x.mean(0) for x in data], 0)
    fig, ax = plt.subplots(figsize=(9, 7.5))         # the layout of 1_window.png (user request 2026-10-07)
    f, u = _unit(times[-1])
    for x in data:
        ax.plot(f * times, np.abs(x.mean(0)), color="gray", lw=0.5, alpha=0.6)
    ax.plot(f * times, np.abs(mu), "k", lw=1.6, label="mean", zorder=3)
    if T2 is not None:
        ax.plot(f * times, t2mod.model(times, T2, n), "C3", lw=2.2, label=f"T2 fit (T2 {f * T2:.3g} {u}, n {n:.2f})")
    ns = data[0].shape[0]
    each = " (ensemble)" if ensemble else f" x {ns} state" + ("s" if ns != 1 else "")
    ax.set(xlabel=f"{TIME_LABEL} ({u})", ylabel=f"|L({TIME_LABEL})|", ylim=L_YLIM, title=f"6. converged coherence: {len(data)} configurations{each}" + (f"\n({fixed})" if fixed else ""))
    ax.grid(alpha=0.3); ax.legend()
    fig.tight_layout(); fig.savefig(out / "6_converged_coherence.png", dpi=120); plt.close(fig)
