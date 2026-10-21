"""Read CCEX coherence curves.

The coherence of one bath state is the `<outfile>_state<s>_wiDiv` file (every inclusion-exclusion exponent applied; the
`_noDiv` file drops the negative exponents and is not the CCE coherence). Each line holds the time and the complex value.

Non-finite values are kept (a numerical failure is reported per configuration, not hidden); CCEX can print them with a space,
e.g. `-inf+inf j`, so the value is read from all tokens after the time.
"""
import hashlib
import re
from pathlib import Path

import numpy as np


def read_curve(path):
    """(times, values) of one CCEX output file. Refuses empty files, short lines, unreadable values and times that are not
    finite or do not strictly increase. Non-finite coherence values are returned as they are (see diagnostics)."""
    t, v = [], []
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < 2:
            raise ValueError(f"{path}:{n}: expected a time and a complex value")
        try:
            t.append(float(parts[0]))
            v.append(complex("".join(parts[1:])))
        except ValueError as e:
            raise ValueError(f"{path}:{n}: {e}") from None
    t, v = np.array(t), np.array(v)
    if len(t) < 2:
        raise ValueError(f"{path}: fewer than 2 time points")
    if not np.all(np.isfinite(t)) or not np.all(np.diff(t) > 0):
        raise ValueError(f"{path}: times must be finite and strictly increasing")
    return t, v


def diagnostics(x, bound=1.01):
    """Per-configuration validity of a state x time array: states with a non-finite value, states with |L| > bound
    (unphysical, U counts both), and the largest finite |L|."""
    finite = np.isfinite(x)
    nonfinite = [int(s) + 1 for s in np.nonzero(~finite.all(1))[0]]
    mag = np.abs(np.where(finite, x, 0))
    over = [int(s) + 1 for s in np.nonzero((mag > bound).any(1))[0]]
    return dict(nonfinite_states=nonfinite, unphysical_states=sorted(set(nonfinite) | set(over)),
                U=len(set(nonfinite) | set(over)), max_abs=float(mag.max()))


def same_grid(a, b):
    """Same time grid up to the output rounding (CCEX prints times with %15g, about 6 significant digits)."""
    if len(a) != len(b):
        return False
    scale = np.maximum(np.abs(a), np.abs(b))
    return bool(np.all(np.abs(a - b) <= 1e-5 * scale + 1e-300))


def read_run(run_dir, prefix="out"):
    """All state curves of one run directory: (times, array state x time). Refuses mixed grids or a missing state."""
    files = list(Path(run_dir).glob(f"{prefix}_state*_wiDiv"))
    if not files:
        raise FileNotFoundError(f"no {prefix}_state*_wiDiv in {run_dir}")
    index = {int(re.search(r"_state(\d+)_wiDiv$", f.name).group(1)): f for f in files}
    ids = sorted(index)
    if ids != list(range(1, len(ids) + 1)):
        raise ValueError(f"{run_dir}: states {ids} are not 1..{len(ids)}")
    times, curves = None, []
    for s in ids:
        t, v = read_curve(index[s])
        if times is None:
            times = t
        elif not same_grid(t, times):
            raise ValueError(f"{index[s]}: time grid differs from state 1")
        curves.append(v)
    return times, np.array(curves)


def _fingerprint(run_dir, prefix):
    h = hashlib.sha256()
    for f in sorted(Path(run_dir).glob(f"{prefix}_state*_wiDiv")):
        h.update(f.read_bytes())
    return h.hexdigest()


def read_runs(run_dirs, prefix="out", expected_states=None, expected_points=None, expected_times=None, manifest_ids=None):
    """One run directory per configuration -> (times, list of state x time arrays, warnings). Configurations may have
    different numbers of states; the time grid must be the same.

    The same directory given twice is refused. Identical output bytes under two names are refused when there is no manifest
    (one configuration would be counted twice and nothing says otherwise); with a manifest (manifest_ids: one unique ID per
    run) they are only a warning, since the manifest, not the bytes, defines the observations. expected_states is an int for
    every run or a list with one count per run (a mixed-size pilot); expected_points / expected_times (deltat * 0..nstep-1)
    must match; a run with fewer states or points, or another grid, is refused (a partial or stale run would otherwise look
    like a smaller complete sample)."""
    warnings = []
    if manifest_ids is not None:
        if len(manifest_ids) != len(run_dirs) or len(set(manifest_ids)) != len(manifest_ids):
            raise ValueError("the manifest needs one unique configuration ID per run")
    seen_path, seen_hash = {}, {}
    for k, d in enumerate(run_dirs):
        key = Path(d).resolve()
        if key in seen_path:
            raise ValueError(f"{d}: the same run directory is given twice ({seen_path[key]})")
        seen_path[key] = d
        fp = _fingerprint(d, prefix)
        if fp in seen_hash:
            if manifest_ids is None:
                raise ValueError(f"{d}: identical outputs to {seen_hash[fp]} (one configuration counted twice?)")
            warnings.append(f"{d} has the same output bytes as {seen_hash[fp]} (manifest IDs differ; check the provenance)")
        seen_hash[fp] = d
    if isinstance(expected_states, (list, tuple)) and len(expected_states) != len(run_dirs):
        raise ValueError("expected_states needs one count per run")
    times, data = None, []
    for k, d in enumerate(run_dirs):
        t, x = read_run(d, prefix)
        want = expected_states[k] if isinstance(expected_states, (list, tuple)) else expected_states
        if want is not None and x.shape[0] != want:
            raise ValueError(f"{d}: {x.shape[0]} states, expected {want} (incomplete run?)")
        if expected_points is not None and len(t) != expected_points:
            raise ValueError(f"{d}: {len(t)} time points, expected {expected_points} (incomplete run?)")
        if expected_times is not None and not same_grid(t, np.asarray(expected_times, float)):
            raise ValueError(f"{d}: time grid differs from the configured one (deltat, nstep)")
        if times is None:
            times = t
        elif not same_grid(t, times):
            raise ValueError(f"{d}: time grid differs from {run_dirs[0]}")
        data.append(x)
    return times, data, warnings
