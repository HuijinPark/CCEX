#!/usr/bin/env python3
"""Read what the C++ writer produced, with python3's json (plan 0B-1, Task 6).

  python json_roundtrip.py test/cluster_verification/check_adaptive/cpp/out/written.json

Every key is checked, by value AND type AND sign, recursively: an earlier draft
checked seven of eleven keys and replacing the rest with "BROKEN" still passed,
and the draft after it still let a nested true become 1, because Python's ==
says True equals 1 (7th and 8th reviews).
"""
import json
import math
from pathlib import Path
import struct
import sys

EXPECTED = {
    # The port's integers are int64 (plan stage 7), so the widest thing the
    # writer can be asked to write is a power of two inside that range - it used
    # to be 2**4096 here, when they were unbounded.
    'work': (int, 2**62),
    'neg': (int, -(2**62)),
    'small': (int, 0),
    'one': (float, 1.),
    'two': (float, 2.),
    'pos_zero': (float, 0.),
    'neg_zero': (float, -0.),
    'tenth': (float, .1),
    'big': (float, 1e308),
    'tiny': (float, 5e-324),
    'inf': (str, 'inf'),
    'ninf': (str, '-inf'),
    'nan': (str, 'nan'),
    'mixed': (list, ['112\\3\n4', True, None]),
    'texts': (dict, {'nul': 'a\x00b', 'accented': '\u00e9', 'emoji': '\U0001f600',
                     'lone': '\ud800', 'literal': '\\u0041'}),
    'dup': (int, 2),
}
ORDER = ['work', 'neg', 'small', 'one', 'two', 'pos_zero', 'neg_zero', 'tenth', 'big', 'tiny',
         'inf', 'ninf', 'nan', 'mixed', 'texts', 'dup']


def same(left, right):
    """Equal in type and value, all the way down.

    type() first: True == 1 in Python, so a nested boolean turned into an integer
    compares equal without it. Floats compare by their bits, which keeps signed
    zero apart and makes NaN equal to NaN, as the comparison rule says.
    """
    if type(left) is not type(right):
        return False
    if isinstance(left, float):
        if left != left and right != right:
            return True                    # NaN equals NaN whatever the payload (design 3.3)
        return struct.pack('<d', left) == struct.pack('<d', right)
    if isinstance(left, list):
        return len(left) == len(right) and all(same(a, b) for a, b in zip(left, right))
    if isinstance(left, dict):
        return list(left) == list(right) and all(same(left[k], right[k]) for k in left)
    return left == right


def main(argv):
    data = json.loads(Path(argv[0]).read_text())
    problems = []
    if list(data) != ORDER:
        problems.append(f'key order {list(data)}')
    for key, (kind, expected) in EXPECTED.items():
        if key not in data:
            problems.append(f'{key} missing')
            continue
        found = data[key]
        if type(found) is not kind:
            problems.append(f'{key} is {type(found).__name__}, expected {kind.__name__}')
        elif not same(found, expected):
            problems.append(f'{key} is {found!r}, expected {expected!r}')
        elif kind is float and math.copysign(1., found) != math.copysign(1., expected):
            problems.append(f'{key} lost its sign')
    for problem in problems:
        print('MISMATCH', problem)
    print(f'{len(EXPECTED)} keys checked')
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
