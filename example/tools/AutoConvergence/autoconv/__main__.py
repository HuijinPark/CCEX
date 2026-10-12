"""Command line: analyse existing CCEX runs (one run directory per bath configuration).

    python -m autoconv analyze RUN_DIR [RUN_DIR ...] [--config DRIVER.json] [--cost-config A --cost-state B] [--prefix out]
                                [--expect-states N] [--expect-points N] [--report-dir DIR]
    python -m autoconv judge --setting RUN_DIR ... --reference RUN_DIR ... [--setting-prefix P] [--reference-prefix P]
                             [--decay-time T] [--name NAME] [--config DRIVER.json] [--report-dir DIR]

analyze: time window and sampling. The sampling step is judged on T2 (2026-10-07): T2 and n from exp(-(t/T2)^n) fitted to
|mean|, its relative SE by linearization (and a two-stage bootstrap), and the cheapest (C, S) for t2_se_target (default 6 %,
the SE scaled by bootstrap/linearized when that is above 1); the curve proxy below stays as a diagnostic.
judge: a truncation setting against a reference (e.g. an rdip_ladder rung, --setting-prefix out_rdip251.359, against the main
output of the same runs), paired by configuration, on [0, T_d] (T_d of the reference unless given, in the output unit):
fit-line difference max <= truncation_max_tol and area <= truncation_area_tol, with the fit check and the smoothed fallback
(autoconv/truncation.py). --report-dir writes PNG figures and text summaries (autoconv/report.py).

Every RUN_DIR is one independent configuration: the same directory or identical outputs twice are refused, and so is a run
with fewer states or time points than expected (from nstate/nstep of base_input with --config, or the options). States
that are non-finite make the report `numerical-failure` (nothing is estimated, no configuration is dropped); finite states
with |L| > 1.01 are listed as warnings. Prints a JSON report: the validity per configuration, the time window (status, decay time, recommended window, per-configuration decay times as
diagnostics), the variance components, sampling_rms_proxy at the present (C, S), the bootstrap norm error
of the present mean, and the cheapest (C, S) for sampling_tolerance when costs are given (bounded by the number of bath files of a fixed
source), and the parts of the total error that are not measured yet. Without --config the defaults of
the auto_convergence specification are used (tolerance 0.01, area norm, decay threshold 0.01, margin 1.5); the sampling error is
judged on [0, T_d] when the time-window step finds a decay time.
"""
import argparse
import json
import sys

import numpy as np

from pathlib import Path

from .config import ConfigError, expected_run_shape, load_config, parse_auto_convergence, read_json, read_manifest, source_capacity
from .curves import diagnostics, read_runs
from .sampling import allocate, bootstrap_norm_error, pooled_within, sampling_rms_proxy, variance_components
from .window import choose_window, decay_time, pooled_mean
from . import report as report_mod
from . import t2 as t2mod
from . import truncation


def _all_finite(*xs):
    return all(np.all(np.isfinite(np.asarray(x, dtype=complex))) for x in xs)


def analyze(run_dirs, settings, cost_config=None, cost_state=None, prefix="out", expected_states=None,
            expected_points=None, max_configurations=None, expected_times=None, manifest_ids=None, max_states=None,
            report_dir=None):
    times, data, warnings = read_runs(run_dirs, prefix, expected_states, expected_points, expected_times, manifest_ids)
    diag = [dict(run=str(d), **diagnostics(x)) for d, x in zip(run_dirs, data)]
    if manifest_ids is not None:
        for g, i in zip(diag, manifest_ids): g["id"] = i
    report = dict(configurations=len(data), states=[int(x.shape[0]) for x in data], time_points=len(times),
                  provenance="manifest" if manifest_ids is not None else "unverified (no manifest)",
                  settings=dict(tolerance=settings.tolerance, error_norm=settings.error_norm,
                                sampling_tolerance=settings.sampling_tolerance,
                                decay_threshold=settings.decay_threshold, window_margin=settings.window_margin),
                  validity=dict(U_total=sum(d["U"] for d in diag), configurations=diag))
    if warnings:
        report["warnings"] = list(warnings)
    failed = [d["run"] for d in diag if d["nonfinite_states"]]
    # Unmeasured parts of the total error are listed, never counted as 0 (plan 0.6 F11, reviews 20261007d R2 and e F9).
    unmeasured = [k for k in ("rdip", "rbath") if not settings.steps[k]]
    if not settings.steps["sampling"]:
        unmeasured.append("sampling (fixed, not assessed)")

    def fail(reason, stage):
        report.update(status="numerical-failure", failed_configurations=failed, reason=reason, failed_stage=stage,
                      total=dict(status="unresolved", unmeasured=unmeasured + ["time_window", "sampling"]))
        return report

    if failed:
        # A numerical failure is not averaged away or dropped: every configuration stays in the report and nothing is
        # estimated from the batch (review 20261007d R8).
        return fail("non-finite states", "input")
    # One sample for every decision (plan 0.8-A, review 20261007e F6): the first S states of each configuration, S = the
    # smallest number of states. Extra states only feed the within-variance check.
    s_common = min(x.shape[0] for x in data)
    balanced = [x[:s_common] for x in data]
    with np.errstate(all="ignore"):
        if settings.steps["time_window"]:
            mu, se = pooled_mean(balanced)
            if not _all_finite(mu, se):
                return fail("non-finite pooled mean or standard error (overflow)", "time_window")
            if settings.window_method == "fit":
                from .window import fit_window
                report["time_window"] = fit_window(times, balanced, settings.decay_threshold, settings.window_margin,
                                                   median_points=settings.window_median_points, points=settings.t2_fit_points)
            else:
                report["time_window"] = choose_window(times, balanced, settings.decay_threshold, settings.window_margin,
                                                      median_points=settings.window_median_points)
        # The area error is judged on the decay part [0, T_d] only (2026-10-07 decision): a long near-zero tail would dilute it.
        mask, evaluated = None, "whole window (no decay time)"
        td = report.get("time_window", {}).get("decay_time")
        if td is not None:
            mask = times <= td * (1 + 1e-9)
            evaluated = f"[0, T_d] = [0, {td:g}] ({int(mask.sum())} of {len(times)} points)"
        if settings.steps["sampling"]:
            comp = variance_components(balanced)
            if not _all_finite(comp["W"], comp["B"], comp["V"]):
                return fail("non-finite variance estimates (overflow)", "sampling")
            c, s = comp["C"], comp["S"]
            sampling = dict(C=c, S=s, states_excluded=int(sum(x.shape[0] for x in data) - c * s), negative_B_points=comp["negative_B"],
                            within_sd_mean=float(np.sqrt(comp["W"]).mean()), between_sd_mean=float(np.sqrt(comp["B"]).mean()),
                            acceptance_statistic="sampling_rms_proxy (an expected-error scale, not a confidence limit)"
                            + ("" if settings.error_norm == "area" else "; max norm: the largest pointwise RMS, not the expected maximum error"),
                            evaluated_on=evaluated,
                            sampling_rms_proxy=sampling_rms_proxy(comp, c, s, settings.error_norm, mask),
                            target=settings.sampling_tolerance)
            if sampling["states_excluded"]:
                w_all = pooled_within(data)
                if not _all_finite(w_all):
                    return fail("non-finite within-configuration variance of all states (overflow)", "sampling")
                sampling["within_sd_mean_all_states"] = float(np.sqrt(w_all).mean())
            boot = bootstrap_norm_error(balanced, norm=settings.error_norm, mask=mask)
            if not _all_finite(boot["mean"], boot["quantile"]):
                return fail("non-finite bootstrap (overflow)", "sampling")
            sampling["bootstrap"] = boot
            sampling["proxy_status"] = "proxy-met" if sampling["sampling_rms_proxy"] <= settings.sampling_tolerance else "sampling-unresolved"
            if cost_config is not None and cost_state is not None:
                sampling["allocation_curve_proxy"] = allocate(comp, settings.sampling_tolerance, cost_config, cost_state,
                                                              norm=settings.error_norm, max_configurations=max_configurations,
                                                              max_states=max_states, mask=mask)
            # The decision (2026-10-07): the relative SE of T2.
            try:
                tc = t2mod.components(times, balanced, points=settings.t2_fit_points)
                boot = t2mod.bootstrap_se_rel(times, balanced, points=settings.t2_fit_points)
            except (ValueError, RuntimeError) as e:          # no decay to fit, or the fit does not converge
                sampling["t2"] = dict(status=f"t2-unresolved: {e}")
                sampling["status"] = "sampling-unresolved"
            else:
                lin = t2mod.se_rel(tc, tc["C"], tc["S"])
                corr = max(1.0, boot / lin) if lin > 0 else 1.0
                t2r = dict(T2=tc["T2"], n=tc["n"], fit_points=tc["fit_points"], B=tc["B"], W=tc["W"], negative_B=bool(tc["negative_B"]),
                           se_rel_linearized=lin, se_rel_bootstrap=boot, target=settings.t2_se_target)
                t2r["status"] = "t2-met" if max(lin, boot) <= settings.t2_se_target else "sampling-unresolved"
                alloc = None
                if cost_config is not None and cost_state is not None:
                    alloc = t2mod.allocate(tc, settings.t2_se_target, cost_config, cost_state, max_configurations=max_configurations,
                                           max_states=max_states, correction=corr)
                    sampling["allocation"] = alloc
                sampling["t2"] = t2r
                sampling["status"] = t2r["status"]
                if report_dir:
                    report_mod.sampling(report_dir, times, balanced, report.get("time_window"), tc, boot, alloc, settings.t2_se_target)
            report["sampling"] = sampling
    if report["validity"]["U_total"]:
        report.setdefault("warnings", []).append(
            f"{report['validity']['U_total']} states with |L| > 1.01 (finite); listed per configuration")
    report["total"] = dict(status="unresolved", unmeasured=unmeasured)
    return report


def _finite(x):
    """JSON-safe copy: numpy scalars to Python numbers, non-finite numbers to null."""
    if isinstance(x, dict):
        return {k: _finite(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_finite(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    return x


# Keys a truncation setting may change between paired runs; every other input key must be equal (review 20261007i).
SETTING_KEYS = {"rdip", "rbath", "rdip_ladder", "outfile"}


def check_pairs(setting_dirs, reference_dirs, unverified=False):
    """Paired runs must be the same configuration: their input.json (CCEX input) must agree on every key but the setting
    keys (the bath, state and axis file paths, seed, physics). Without an input.json the pairing is refused unless
    unverified is set (then it is reported as a warning)."""
    if len(setting_dirs) != len(reference_dirs):
        raise ValueError("--setting and --reference need the same number of runs")
    warnings = []
    for s, r in zip(setting_dirs, reference_dirs):
        fs, fr = Path(s) / "input.json", Path(r) / "input.json"
        if not (fs.is_file() and fr.is_file()):
            if not unverified:
                raise ValueError(f"{s} / {r}: no input.json to check that they are the same configuration "
                                 "(--unverified-pairs to accept the order as the pairing)")
            warnings.append(f"{s} / {r}: pairing not verified (no input.json)")
            continue
        a, b = read_json(fs), read_json(fr)
        diff = sorted(k for k in set(a) | set(b) if k not in SETTING_KEYS and a.get(k) != b.get(k))
        if diff:
            raise ValueError(f"{s} / {r}: not the same configuration, inputs differ in {diff}")
    return warnings


def judge_runs(setting_dirs, reference_dirs, settings, setting_prefix="out", reference_prefix="out", decay=None,
               name="setting", report_dir=None, unverified_pairs=False):
    """Truncation judgement of paired runs (configuration k of the setting against configuration k of the reference)."""
    pair_warnings = check_pairs(setting_dirs, reference_dirs, unverified_pairs)
    ts, ds, _ = read_runs(setting_dirs, setting_prefix)
    tr, dr, _ = read_runs(reference_dirs, reference_prefix)
    if len(ts) != len(tr) or not np.allclose(ts, tr, rtol=1e-5, atol=0):
        raise ValueError("setting and reference have different time grids")
    for a, b, d in zip(ds, dr, setting_dirs):
        if a.shape != b.shape:
            raise ValueError(f"{d}: {a.shape[0]} states in the setting, {b.shape[0]} in the reference")
    if not all(_all_finite(x) for x in ds + dr):
        return dict(verdict="numerical-failure", reason="non-finite state values")
    if decay is None:
        mu, se = pooled_mean(dr)
        decay, _ = decay_time(tr, mu, se, settings.decay_threshold, median_points=settings.window_median_points)
        if decay is None:
            return dict(verdict="window-unresolved", reason="the reference does not decay inside its window; give --decay-time")
    j = truncation.judge(tr, np.array([x.mean(0) for x in ds]), np.array([x.mean(0) for x in dr]), decay,
                         max_tol=settings.truncation_max_tol, area_tol=settings.truncation_area_tol, points=settings.t2_fit_points,
                         sigma=decay / 25)
    if pair_warnings:
        j["warnings"] = pair_warnings
    if report_dir:
        same = setting_prefix == reference_prefix
        report_mod.truncation(report_dir, name, tr, j, "setting" if same else setting_prefix, "reference" if same else reference_prefix)
    j.pop("_curves")
    return j


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == "run":                      # main.out --auto_conv hands over here (autoconv/runner.py)
        from . import runner
        return runner.main(argv[1:])
    if argv and argv[0] == "status":                   # main.out --auto_status (where a run is, without computing)
        from . import runner
        return runner.status_main(argv[1:])
    ap = argparse.ArgumentParser(prog="autoconv", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    a = sub.add_parser("analyze")
    a.add_argument("runs", nargs="*")
    a.add_argument("--manifest", help="JSON manifest: configuration IDs, run folders, state counts, available configurations, grid")
    a.add_argument("--config")
    a.add_argument("--cost-config", type=float)
    a.add_argument("--cost-state", type=float)
    a.add_argument("--prefix", default="out")
    a.add_argument("--expect-states", type=int, help="states every run must have (default: nstate of base_input)")
    a.add_argument("--expect-points", type=int, help="time points every run must have (default: nstep of base_input)")
    a.add_argument("--report-dir", help="write sampling.png and sampling.txt here")
    j = sub.add_parser("judge")
    j.add_argument("--setting", nargs="+", required=True)
    j.add_argument("--reference", nargs="+", required=True)
    j.add_argument("--setting-prefix", default="out")
    j.add_argument("--reference-prefix", default="out")
    j.add_argument("--decay-time", type=float, help="T_d in the output unit (default: from the reference)")
    j.add_argument("--name", default="setting")
    j.add_argument("--config")
    j.add_argument("--report-dir", help="write truncation_<name>.png and .txt here")
    j.add_argument("--unverified-pairs", action="store_true", help="accept runs without input.json, paired by order (warned)")
    args = ap.parse_args(argv)
    if args.command == "judge":
        try:
            settings = load_config(args.config) if args.config else parse_auto_convergence({})
            rep = judge_runs(args.setting, args.reference, settings, args.setting_prefix, args.reference_prefix,
                             args.decay_time, args.name, args.report_dir, args.unverified_pairs)
        except (ConfigError, ValueError, FileNotFoundError, RuntimeError) as e:
            print(f"autoconv: {e}", file=sys.stderr)
            return 2
        json.dump(_finite(rep), sys.stdout, indent=1, allow_nan=False)
        print()
        return 0
    try:
        settings = load_config(args.config) if args.config else parse_auto_convergence({})
        if (args.cost_config is None) != (args.cost_state is None):
            raise ConfigError("give both --cost-config and --cost-state, or neither")
        states, points = expected_run_shape(settings) if args.config else (None, None)
        times, ids, capacity, max_states = None, None, source_capacity(settings) if args.config else None, None
        runs = args.runs
        if args.manifest:
            if runs:
                raise ConfigError("give the runs either on the command line or in the manifest, not both")
            runs, ids, states, available, grid, max_states = read_manifest(args.manifest)
            if available is not None:
                capacity = int(available)
            if grid:
                points = int(grid["nstep"])
                times = float(grid["deltat"]) * np.arange(points)
        elif not runs:
            raise ConfigError("no runs given")
        if times is None and args.config and points and "deltat" in settings.base_input:
            times = float(settings.base_input["deltat"]) * np.arange(points)
        report = analyze(runs, settings, args.cost_config, args.cost_state, args.prefix,
                         args.expect_states if args.expect_states is not None else states,
                         args.expect_points if args.expect_points is not None else points,
                         capacity, times, ids, max_states, args.report_dir)
    except (ConfigError, ValueError, FileNotFoundError, OverflowError) as e:
        print(f"autoconv: {e}", file=sys.stderr)
        return 2
    json.dump(_finite(report), sys.stdout, indent=1, allow_nan=False)
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
