"""Compare two refinement JSON files after removing elapsed measurements.

Usage: python compare_summaries.py A/summary.json B/summary.json
Only keys that hold an elapsed measurement are dropped (`wall_s`, `search_s`,
`seconds_used`, `timing` and any `*_wall_s`), and only when asked. Limits,
settings and decision flags such as `max_total_seconds`, `seconds_limit` or
`raw_propagation_started` are compared, and requests, inputs and native
responses are compared with nothing dropped. Each run's own directory is
replaced by <ROOT>. Prints the differing paths; exit status 1 if any, 2 if a
file is missing or is not valid JSON. A comparison that did not happen is never
"no difference".
"""
import json
import math
from pathlib import Path
import sys

MEASURED = frozenset({'wall_s', 'search_s', 'seconds_used', 'timing'})


class LoadError(Exception):
    """A file could not be read as JSON."""


def measured(key):
    """True only for keys holding an elapsed measurement, never a setting."""
    return key in MEASURED or key.endswith('_wall_s')


def normalize(x, root, drop=False):
    if isinstance(x, dict):
        return {k: normalize(v, root, drop) for k, v in x.items() if not (drop and measured(k))}
    if isinstance(x, list):
        return [normalize(v, root, drop) for v in x]
    return x.replace(root, '<ROOT>') if isinstance(x, str) else x


def load(path, root=None, drop=False):
    """Normalized JSON content; raises LoadError instead of returning nothing."""
    path = Path(path)
    try:
        data = json.loads(path.read_text())
    except (OSError, ValueError) as error:
        raise LoadError(f'{path}: {error}') from error
    return normalize(data, str(path.parent if root is None else root), drop)


def same(a, b):
    """Python equality, except that -0.0 differs from 0.0 and NaN equals NaN.

    A report writes non-finite values as strings, so a bare NaN here means two
    runs produced the same undefined value, not a difference. A sign flip on
    zero is a real difference and the design requires reporting it.
    """
    if isinstance(a, float) and isinstance(b, float):
        if math.isnan(a) or math.isnan(b):
            return math.isnan(a) and math.isnan(b)
        return a == b and math.copysign(1., a) == math.copysign(1., b)
    return a == b


def diff(a, b, path=''):
    if type(a) is not type(b):
        yield path
    elif isinstance(a, dict):
        for k in sorted(a.keys() | b.keys()):
            yield from ([f'{path}/{k}'] if k not in a or k not in b else diff(a[k], b[k], f'{path}/{k}'))
    elif isinstance(a, list):
        if len(a) != len(b):
            yield path
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                yield from diff(x, y, f'{path}[{i}]')
    elif not same(a, b):
        yield path


def compare(a, b, root=None, drop=False):
    return list(diff(load(a, root, drop), load(b, root, drop)))


def main(argv):
    try:
        bad = compare(argv[0], argv[1], drop=True)
    except LoadError as error:
        print(error, file=sys.stderr)
        return 2
    print('\n'.join(bad) or 'identical after normalization')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
