"""What the analyzer knows about a CCEX coherence file, without the plotting.

This module needs numpy and nothing else, so it can be imported from any script:

    import sys; sys.path.insert(0, 'example/analysis/Analyzer/CoherenceAnalyzer.v1')
    from coherence import read_coherence, error_correct, ensemble_average, write_coherence

    t, L_wD = read_coherence('CCE2_..._conf1_wiDiv')
    _, L_nD = read_coherence('CCE2_..._conf1_noDiv')
    L_EC = error_correct(L_wD, L_nD, tolerance=2.0)
    write_coherence('CCE2_..._conf1_ErrCorr', t, L_EC)

The file format: one line per time step, the time in ms then the complex
coherence as "+re+imj" (or as two real columns), tab-separated.

The error-correction rule (the group's convention, not a physical law; the
regression references were made with it, so do not change it):
  CCE divides the coherence of each cluster by that of its sub-clusters, and the
  division blows up where a sub-cluster coherence passes through zero. The divided
  curve (wiDiv) is therefore checked against the undivided one (noDiv). A point
  where |wiDiv - noDiv| >= tolerance is a divergence: the first and last points
  become 0, any other point becomes the mean of the previous corrected value and
  the next non-divergent wiDiv value. If the last point diverges it is first pulled
  to noDiv + tolerance - 1e-6 so the search for a "next good point" terminates.
`DataSpace.py` and `run.py` call these functions; they are the one implementation.
"""
import itertools
import os

import numpy as np

SPLITER = " | "


def _rows(path):
    """Whitespace-split tokens per line. A missing file is reported and gives no rows."""
    rows = []
    try:
        with open(path, 'r') as f:
            for line in f:
                rows.append(line.strip().split())
    except FileNotFoundError:
        print(f"Error ! {path} file not found !")
    return rows


def read_coherence(path, time_scale=1.0, stride=1):
    """Read one coherence file: (t, L) with t float64 (multiplied by time_scale) and
    L complex. Accepts "t  re+imj" and "t  re  im"; a "nan" in the coherence column
    that keeps the file from parsing is read as 0. Raises ValueError if the file
    cannot be read at all."""
    rows = _rows(path)
    try:
        data = np.array(rows).astype(complex)
    except (ValueError, TypeError):
        for row in rows:
            if len(row) > 1 and 'nan' in row[1]:
                row[:] = [row[0], 0.0]
        data = np.array(rows).astype(complex)
    if data.ndim != 2 or data.shape[1] not in (2, 3):
        raise ValueError(f"{path}: expected 2 or 3 columns")
    t = data[:, 0].real * time_scale
    L = data[:, 1] + data[:, 2] * 1j if data.shape[1] == 3 else data[:, 1]
    return t[::stride], L[::stride]


def writefile(fout, dtype, data, comment=None):
    """Columns to a text file. data = [col, col, ...], dtype one letter per column:
    f  fixed 10 decimals   d  integer   s  string   c  complex as +re+imj   e  exponent."""
    with open(fout, 'w') as f:
        for col in range(len(data[0])):
            for row in range(len(data)):
                if dtype[row] == 'f':
                    line = "{:>.10f}\t".format(data[row][col])
                elif dtype[row] == 'd':
                    line = "{:>10d}\t".format(data[row][col])
                elif dtype[row] == 's':
                    line = "{:>10s}\t".format(data[row][col])
                elif dtype[row] == 'c':
                    line = "{:>+.10f}{:>+.10f}j\t".format(data[row][col].real, data[row][col].imag)
                elif dtype[row] == 'e':
                    line = "{:>.10E}\t".format(data[row][col])
                else:
                    raise ValueError(f"writefile: unknown column type {dtype[row]!r}")
                f.write(line)
            f.write("\n")
        if comment is not None:
            f.write(comment)


def write_coherence(path, t, L):
    """The file format CCEX and the analyzer share: time, then +re+imj."""
    writefile(path, ['f', 'c'], [np.asarray(t).real, np.asarray(L, complex)])


def error_correct(L_wD, L_nD, tolerance=2.0):
    """The rule in the module docstring. L_wD is the divided curve, L_nD the undivided
    reference; both complex arrays over the same time grid. Like the original, this
    may pull the last element of L_wD in place, and it prints each corrected point."""
    L_EC = []
    n_step = len(L_wD)
    if n_step == 0:
        return np.asarray(L_EC, complex)

    # ">=" as in the rule above and in the loop below: with ">" a last point exactly at the
    # tolerance was not pulled and the search for a next good point ran past the end (D25)
    if (np.abs(L_wD[n_step-1] - L_nD[n_step-1]) >= tolerance):
        L_wD[n_step-1] = L_nD[n_step-1] + tolerance - 0.000001

    for i in range(len(L_wD)):
        error = np.abs(L_wD[i] - L_nD[i])

        if error >= tolerance:
            if i == 0 or i == (n_step - 1):
                L_EC.append(0)
            else:
                count = 1
                is_end = True
                while (is_end):
                    if (np.abs(L_wD[i+count] - L_nD[i+count]) < tolerance):
                        L_EC.append((L_EC[i-1] + L_wD[i+count])/2)
                        print(f"\t\tDo error correction for {i:5}th coherence : {L_wD[i]:30.5f}   ->   {L_EC[-1]:15.5f}")
                        is_end = False
                    count = count+1
        else:
            L_EC.append(L_wD[i])

    return np.asarray(L_EC, complex)


def ensemble_average(curves):
    """Mean of several (t, L) pairs over the same time grid. Returns (t, mean, count).
    Nothing passed in is modified."""
    curves = list(curves)
    if not curves:
        raise ValueError("ensemble_average: no curves")
    t0 = np.asarray(curves[0][0])
    total = np.zeros(len(t0), complex)
    for t, L in curves:
        if not np.array_equal(np.asarray(t), t0):
            raise ValueError("ensemble_average: the time grids differ")
        total += np.asarray(L, complex)
    return t0, total / len(curves), len(curves)


def expand_names(directory, pattern, variables):
    """The analyzer's file-naming convention. variables is [(NAME, [value, ...]), ...];
    every NAME in pattern is replaced by each of its values, in order, and the key of a
    file is its values joined by " | ". Returns {key: path}."""
    names = [name for name, _ in variables]
    out = {}
    for values in itertools.product(*[list(map(str, v)) for _, v in variables]):
        path = os.path.join(directory, pattern) if directory else pattern
        for name, value in zip(names, values):
            path = path.replace(name, value)
        out[SPLITER.join(values)] = path
    return out
