#!/usr/bin/env python3
"""
CCEX regression comparison
==========================

Compares the coherences produced by a modified build against a reference
version and draws the RMS error as a grouped bar chart.

A "version folder" is any subdirectory of test/regression/ that holds the
per-test layout produced by test.sh:

    <version>/
        test.sh
        <test>/job
        <test>/inputfiles/ccein_*.json
        <test>/results/<nstate>/raw/<coherence files>

v1.0.0 is the reference. To check a code change, in one command:

    ./compare.py -n +qgyroFIX -x   # clone, run, compare
                                   # -> rms_+qgyroFIX_vs_v1.0.0.png

or step by step, which is the same thing:

    ./compare.py -n +qgyroFIX      # clone the v1.0.0 layout, no results
    cd +qgyroFIX && ./test.sh      # run it against the modified build
    cd .. && ./compare.py          # -> rms_+qgyroFIX_vs_v1.0.0.png

-x never touches the reference: running it on the reference version is
refused, so a stray ./compare.py v1.0.0 -x cannot overwrite the frozen
results everything else is measured against.

Release tags (v1.0.1, v1.1.0, ...) are reserved for what lands on main.
While a fix or an optimization is still on a branch, name its folder after
the change itself with a leading '+', read as "reference plus this":

    +qgyroFIX          v1.0.0 with the qgyro float->double fix
    +eigPropagator     v1.0.0 with the eigendecomposition propagator

Any folder name works; the '+' form is just the convention.

Usage
-----
    ./compare.py                       compare every non-reference version
    ./compare.py <target> [<target>..] compare specific versions
    ./compare.py -n <name>             clone the reference layout (no results)
    ./compare.py -x                    run each target's test.sh first
    ./compare.py -r <ref>              use another reference (default v1.0.0)
    ./compare.py -l                    list version folders
    ./compare.py -h                    this text

Output
------
Two figures per comparison, plus the tables they were drawn from:

    rms_<target>_vs_<ref>.png     accuracy -- did the numbers move?
    wall_<target>_vs_<ref>.png    cost     -- did it get faster?

Accuracy is one RMS per (calculation, quantity), always over the COMPLEX
coherence:

    RMS = sqrt( sum_i |L_ref(t_i) - L_new(t_i)|^2 / N )

A calculation is one test at one nstate.  The four quantities are the raw
CCEX outputs _noDiv and _wiDiv, the analyzer's error-corrected _ErrCorr,
and the ensemble average EnsAvg.  When a calculation has several bath
configurations, all of their points are pooled into that one RMS.

Cost is the wall time from the '[run] ... wall=..s' footer each job writes
into its .process file, summed over the configurations of a calculation.
Those files can grow to many GB, so only their last 32 kB is ever read.
The figure shows absolute times on a log axis and the speedup ref/new per
calculation; it is one run per calculation, so treat small ratios as noise.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = Path(__file__).resolve().parent
REF_DEFAULT = "v1.0.0"

# order matters: it is the bar order inside a group
QUANTITIES = ["noDiv", "wiDiv", "ErrCorr", "EnsAvg"]
COLORS = {
    "noDiv": "#4C72B0",
    "wiDiv": "#DD8452",
    "ErrCorr": "#55A868",
    "EnsAvg": "#C44E52",
}
# bars this size or smaller are treated as "no difference at all"
EXACT_ZERO = 0.0


# --------------------------------------------------------------------------
# discovery
# --------------------------------------------------------------------------
def version_dirs():
    """Subdirectories that look like a regression version folder."""
    out = []
    for d in sorted(p for p in HERE.iterdir() if p.is_dir()):
        if any((t / "job").is_file() for t in d.iterdir() if t.is_dir()):
            out.append(d.name)
    return out


def calculations(version):
    """[(test, nstate, raw_dir), ...] for every calculation that has results."""
    root = HERE / version
    found = []
    for test in sorted(p for p in root.iterdir() if p.is_dir()):
        results = test / "results"
        if not results.is_dir():
            continue
        for nst in sorted(p for p in results.iterdir() if p.is_dir()):
            raw = nst / "raw"
            if raw.is_dir():
                found.append((test.name, nst.name, raw))
    return found


def classify(name):
    if name.endswith("_noDiv"):
        return "noDiv"
    if name.endswith("_wiDiv"):
        return "wiDiv"
    if name.endswith("_ErrCorr"):
        return "ErrCorr"
    if "EnsAvg" in name:
        return "EnsAvg"
    return None


# --------------------------------------------------------------------------
# numerics
# --------------------------------------------------------------------------
def load_coherence(path):
    """Read '<time>\\t<+re+imj>' rows. Returns (t, L) as float / complex arrays."""
    t, z = [], []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                tt = float(parts[0])
                zz = complex(parts[1])
            except ValueError:
                continue  # header or comment
            t.append(tt)
            z.append(zz)
    return np.asarray(t), np.asarray(z, dtype=complex)


RUN_FOOTER = re.compile(r"\[run\]\s+exit=(\d+)\s+cores=(\d+)\s+wall=([0-9.]+)s")


def read_wall(path, tail_bytes=32768):
    """Last '[run] exit=.. cores=.. wall=..s' footer of a .process file.

    Seek-based on purpose: a .process file can reach many GB, and the footer is
    always at the end, so never read the whole thing.
    """
    size = os.path.getsize(path)
    with open(path, "rb") as fh:
        fh.seek(max(0, size - tail_bytes))
        tail = fh.read().decode("utf-8", "replace")
    last = None
    for m in RUN_FOOTER.finditer(tail):
        last = m
    if last is None:
        return None
    return {"exit": int(last.group(1)), "cores": int(last.group(2)),
            "wall": float(last.group(3))}


def wall_for_calculation(nstate_dir):
    """Total wall over every .process of one calculation, plus the core counts."""
    total, cores, n = 0.0, set(), 0
    for p in sorted(nstate_dir.glob("*.process")):
        got = read_wall(p)
        if got is None:
            continue
        total += got["wall"]
        cores.add(got["cores"])
        n += 1
    return (total, cores, n) if n else (None, cores, 0)


def rms_for_calculation(ref_raw, tgt_raw):
    """{quantity: (rms, npoints, note)} pooling every file of that quantity."""
    acc = {q: {"sq": 0.0, "n": 0, "note": ""} for q in QUANTITIES}

    for ref_file in sorted(ref_raw.iterdir()):
        if not ref_file.is_file():
            continue
        q = classify(ref_file.name)
        if q is None:
            continue

        tgt_file = tgt_raw / ref_file.name
        if not tgt_file.is_file():
            acc[q]["note"] = "missing"
            continue

        _, a = load_coherence(ref_file)
        _, b = load_coherence(tgt_file)
        if a.shape != b.shape:
            acc[q]["note"] = f"shape {a.shape[0]} vs {b.shape[0]}"
            continue

        d = a - b
        acc[q]["sq"] += float(np.sum(np.abs(d) ** 2))
        acc[q]["n"] += d.size

    out = {}
    for q, v in acc.items():
        if v["n"] == 0:
            out[q] = (None, 0, v["note"] or "no file")
        else:
            out[q] = (float(np.sqrt(v["sq"] / v["n"])), v["n"], v["note"])
    return out


def compare(target, ref):
    """[(test, nstate, {quantity: (rms, n, note)}, wall_ref, wall_tgt, cores), ...]"""
    ref_calcs = calculations(ref)
    if not ref_calcs:
        sys.exit(f"ERROR: reference '{ref}' has no results. Run its test.sh first.")

    rows = []
    for test, nstate, ref_raw in ref_calcs:
        tgt_raw = HERE / target / test / "results" / nstate / "raw"
        wr, cr, _ = wall_for_calculation(ref_raw.parent)
        if not tgt_raw.is_dir():
            rows.append((test, nstate,
                         {q: (None, 0, "not run") for q in QUANTITIES}, wr, None, cr))
            continue
        wt, ct, _ = wall_for_calculation(tgt_raw.parent)
        rows.append((test, nstate, rms_for_calculation(ref_raw, tgt_raw),
                     wr, wt, cr | ct))
    return rows


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------
def calc_label(test, nstate):
    return f"{test}/{nstate}"


def print_table(rows, target, ref):
    wide = max(len(calc_label(t, n)) for t, n, *_ in rows)
    head = f"  {'calculation':<{wide}}  " + "  ".join(f"{q:>11}" for q in QUANTITIES)
    print(f"\nRMS |dL|   {target} vs {ref}\n")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for test, nstate, res, *_ in rows:
        cells = []
        for q in QUANTITIES:
            rms, _, note = res[q]
            if rms is None:
                cells.append(f"{note:>11}")
            elif rms == EXACT_ZERO:
                cells.append(f"{'0 (exact)':>11}")
            else:
                cells.append(f"{rms:>11.3e}")
        print(f"  {calc_label(test, nstate):<{wide}}  " + "  ".join(cells))

    vals = [r[0] for _, _, res, *_ in rows for r in res.values() if r[0] is not None]
    if vals:
        print(f"\n  max RMS = {max(vals):.3e}   ({sum(v == 0 for v in vals)}"
              f"/{len(vals)} comparisons are exactly 0)")


def print_wall_table(rows, target, ref):
    usable = [r for r in rows if r[3] and r[4]]
    if not usable:
        print("\n  (no wall times: .process footers missing on one side)")
        return None

    wide = max(len(calc_label(t, n)) for t, n, *_ in rows)
    print(f"\nwall time [s]   {target} vs {ref}\n")
    print(f"  {'calculation':<{wide}}  {ref[:11]:>11}  {target[:11]:>11}  "
          f"{'speedup':>8}")
    print("  " + "-" * (wide + 36))
    for test, nstate, _, wr, wt, cores in rows:
        lab = calc_label(test, nstate)
        if not wr or not wt:
            print(f"  {lab:<{wide}}  {'-' if not wr else f'{wr:.2f}':>11}  "
                  f"{'-' if not wt else f'{wt:.2f}':>11}  {'-':>8}")
            continue
        print(f"  {lab:<{wide}}  {wr:>11.2f}  {wt:>11.2f}  {wr / wt:>7.3f}x")

    tr = sum(r[3] for r in usable)
    tt = sum(r[4] for r in usable)
    print("  " + "-" * (wide + 36))
    print(f"  {'total':<{wide}}  {tr:>11.2f}  {tt:>11.2f}  {tr / tt:>7.3f}x")

    allcores = set().union(*(r[5] for r in usable))
    if len(allcores) > 1:
        print(f"\n  !! core counts differ across runs ({sorted(allcores)})"
              f" -- the comparison is not apples to apples")
    print("  note: single run per calculation, so small ratios are within noise")
    return tr / tt


def short_labels(rows):
    return [f"{t.replace('Diamond_NV_', '')}\n{n}" for t, n, *_ in rows]


def plot(rows, target, ref, out_png):
    labels = short_labels(rows)
    ngroup, nbar = len(rows), len(QUANTITIES)
    x = np.arange(ngroup)
    width = 0.8 / nbar

    vals = [r[0] for _, _, res, *_ in rows for r in res.values() if r[0] is not None]
    nonzero = [v for v in vals if v > 0]
    if nonzero:
        floor = 10 ** (np.floor(np.log10(min(nonzero))) - 1)
        ceil = 10 ** (np.ceil(np.log10(max(nonzero))) + 1)
    else:
        floor, ceil = 1e-17, 1e-12

    fig, ax = plt.subplots(figsize=(1.6 * ngroup + 3, 6.0))

    for k, q in enumerate(QUANTITIES):
        pos = x - 0.4 + width * (k + 0.5)
        for i, (_, _, res, *_rest) in enumerate(rows):
            rms, _, note = res[q]
            if rms is None:
                ax.text(pos[i], floor * 2, note, rotation=90, ha="center",
                        va="bottom", fontsize=6, color="0.45")
                continue
            if rms == EXACT_ZERO:
                # log axis cannot show 0 -- mark it on the floor instead
                ax.plot(pos[i], floor * 1.6, marker="D", ms=5, mfc="none",
                        mec=COLORS[q], mew=1.4, clip_on=False, zorder=5)
                continue
            ax.bar(pos[i], rms, width * 0.92, bottom=floor,
                   color=COLORS[q], zorder=3)
            ax.text(pos[i], rms * 1.25, f"{rms:.1e}", rotation=90, ha="center",
                    va="bottom", fontsize=6.5, color="0.25")

    handles = [Patch(facecolor=COLORS[q], label=q) for q in QUANTITIES]
    handles.append(Line2D([], [], marker="D", ms=5, ls="none", mfc="none",
                          mec="0.35", label="0 (exact)"))

    if vals and not nonzero:
        # every comparison is bit-identical -- say so, the axis alone looks empty
        ax.text(0.5, 0.5, f"no difference\nall {len(vals)} comparisons are exactly 0",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=15, color="#55A868", weight="bold")

    ax.set_yscale("log")
    ax.set_ylim(floor, ceil)
    ax.set_xlim(-0.6, ngroup - 0.4)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel(r"RMS $|\Delta L|$   $\sqrt{\sum_i |L_{ref}-L_{new}|^2/N}$")
    ax.set_title(f"CCEX regression : complex coherence RMS error\n"
                 f"{target}  vs  {ref}", fontsize=11)
    ax.grid(axis="y", which="both", ls=":", lw=0.6, color="0.75", zorder=0)
    ax.set_axisbelow(True)
    ax.legend(handles=handles, ncol=len(handles), fontsize=8,
              loc="upper center", bbox_to_anchor=(0.5, -0.13), frameon=False)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_png, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_wall(rows, target, ref, out_png):
    """Two panels: absolute wall time side by side, and speedup per calculation."""
    usable = [r for r in rows if r[3] and r[4]]
    if not usable:
        return False

    labels = short_labels(usable)
    x = np.arange(len(usable))
    wr = np.array([r[3] for r in usable])
    wt = np.array([r[4] for r in usable])
    speed = wr / wt

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(1.6 * len(usable) + 3, 8.0),
        gridspec_kw={"height_ratios": [2, 1.4], "hspace": 0.35})

    # --- absolute wall time ---
    w = 0.38
    ax1.bar(x - w / 2, wr, w, label=ref, color="#8C8C8C", zorder=3)
    ax1.bar(x + w / 2, wt, w, label=target, color="#4C72B0", zorder=3)
    for xi, (a, b) in enumerate(zip(wr, wt)):
        ax1.text(xi - w / 2, a, f"{a:.2f}", ha="center", va="bottom", fontsize=7)
        ax1.text(xi + w / 2, b, f"{b:.2f}", ha="center", va="bottom", fontsize=7)
    ax1.set_yscale("log")
    ax1.set_ylim(min(wr.min(), wt.min()) / 3, max(wr.max(), wt.max()) * 4)
    ax1.set_ylabel("wall time [s]  (log)")
    ax1.set_xticks(x)
    ax1.set_xticklabels([""] * len(x))
    ax1.legend(fontsize=8, ncol=2, frameon=False, loc="upper left")
    ax1.grid(axis="y", which="both", ls=":", lw=0.6, color="0.75", zorder=0)
    ax1.set_axisbelow(True)

    tr, tt = wr.sum(), wt.sum()
    ax1.set_title(f"CCEX regression : wall time\n{target}  vs  {ref}"
                  f"      total {tr:.2f}s -> {tt:.2f}s  ({tr / tt:.3f}x)",
                  fontsize=11)

    # --- speedup ---
    colors = ["#55A868" if s > 1.02 else "#C44E52" if s < 0.98 else "#8C8C8C"
              for s in speed]
    ax2.bar(x, speed - 1.0, 0.6, bottom=1.0, color=colors, zorder=3)
    for xi, s in enumerate(speed):
        ax2.text(xi, s, f"{s:.3f}x", ha="center",
                 va="bottom" if s >= 1 else "top", fontsize=7.5)
    ax2.axhline(1.0, color="0.3", lw=1.2, zorder=4)
    span = max(0.05, float(np.abs(speed - 1.0).max()) * 1.6)
    ax2.set_ylim(1 - span, 1 + span)
    ax2.set_ylabel("speedup  (ref / new)")
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, fontsize=8)
    ax2.grid(axis="y", ls=":", lw=0.6, color="0.75", zorder=0)
    ax2.set_axisbelow(True)
    ax2.text(0.995, 0.04, "above 1.0 = faster   |   single run, small "
             "deviations are noise", transform=ax2.transAxes, ha="right",
             va="bottom", fontsize=7.5, color="0.45")

    for ax in (ax1, ax2):
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    fig.savefig(out_png, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return True


# --------------------------------------------------------------------------
# new-version scaffolding
# --------------------------------------------------------------------------
SKIP_NAMES = {"results", "job.log"}


def clone_version(name, source):
    src, dst = HERE / source, HERE / name
    if not src.is_dir():
        sys.exit(f"ERROR: source version '{source}' not found in {HERE}")
    if dst.exists():
        sys.exit(f"ERROR: '{name}' already exists. Remove it first.")

    dst.mkdir()
    made = []

    top = src / "test.sh"
    if top.is_file():
        shutil.copy2(top, dst / "test.sh")
        made.append("test.sh")

    for test in sorted(p for p in src.iterdir() if p.is_dir()):
        if not (test / "job").is_file():
            continue
        (dst / test.name).mkdir()
        for item in sorted(test.iterdir()):
            if item.name in SKIP_NAMES or item.name.endswith(".swp"):
                continue
            if item.is_dir():
                shutil.copytree(item, dst / test.name / item.name)
            else:
                shutil.copy2(item, dst / test.name / item.name)
            made.append(f"{test.name}/{item.name}")

    print(f"created {name}/ from {source}  (results/ and logs excluded)")
    for m in made:
        print(f"    {m}")


def run_tests(version, ref):
    """Run <version>/test.sh. Refuses to touch the reference."""
    if version == ref:
        sys.exit(f"ERROR: refusing to run '{version}' -- it is the reference.\n"
                 f"       Re-running it would overwrite the frozen results that\n"
                 f"       every other version is measured against.")
    script = HERE / version / "test.sh"
    if not script.is_file():
        sys.exit(f"ERROR: {version}/test.sh not found")

    print(f"\n{'=' * 60}\n>> running {version}/test.sh\n{'=' * 60}", flush=True)
    rc = subprocess.call([str(script)], cwd=str(HERE / version))
    if rc != 0:
        print(f"\n!! {version}/test.sh exited {rc} -- comparing anyway",
              file=sys.stderr)
    return rc


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        prog="compare.py", add_help=False,
        description="Compare CCEX regression coherences against a reference.")
    ap.add_argument("targets", nargs="*")
    ap.add_argument("-n", "--new", metavar="NAME")
    ap.add_argument("-x", "--run", action="store_true")
    ap.add_argument("-r", "--ref", default=REF_DEFAULT, metavar="REF")
    ap.add_argument("-l", "--list", action="store_true")
    ap.add_argument("-h", "--help", action="store_true")
    ns = ap.parse_args()

    if ns.help:
        print(__doc__)
        return 0

    versions = version_dirs()

    if ns.list:
        print(f"version folders in {HERE.name}/:")
        wide = max((len(v) for v in versions), default=1)
        for v in versions:
            n = len(calculations(v))
            tag = "  (reference)" if v == ns.ref else ""
            print(f"    {v:<{wide}}   {n} calculation(s) with results{tag}")
        return 0

    ref = ns.ref

    if ns.new:
        clone_version(ns.new, ref)
        if not ns.run:
            print(f"\nnext:  ./'{ns.new}'/test.sh && ./compare.py"
                  f"      (or re-run this with -x)")
            return 0
        targets = [ns.new]
        versions = version_dirs()
    elif ns.targets:
        targets = ns.targets
    else:
        targets = [v for v in versions if v != ref]
        if not targets:
            sys.exit(f"ERROR: no version to compare. Reference is '{ref}'.\n"
                     f"       Make one with:  ./compare.py -n <name> -x")

    rc = 0
    if ns.run:
        for target in targets:
            if run_tests(target, ref) != 0:
                rc = 1
    for target in targets:
        if target not in versions:
            print(f"ERROR: '{target}' is not a version folder. "
                  f"available: {', '.join(versions)}", file=sys.stderr)
            rc = 1
            continue
        rows = compare(target, ref)
        root = HERE.parent.parent.parent

        print_table(rows, target, ref)
        out = HERE / f"rms_{target}_vs_{ref}.png"
        plot(rows, target, ref, out)
        print(f"\n  wrote {out.relative_to(root)}")

        print_wall_table(rows, target, ref)
        outw = HERE / f"wall_{target}_vs_{ref}.png"
        if plot_wall(rows, target, ref, outw):
            print(f"\n  wrote {outw.relative_to(root)}\n")
        else:
            print()
    return rc


if __name__ == "__main__":
    sys.exit(main())
