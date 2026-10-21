#!/usr/bin/env python3
"""Compare the tree C++ read with the tree CPython read (plan 0B-1, Task 7).

  python json_expect.py cases  <out-dir> <json-cases-dir>
  python json_expect.py tree   <out-dir> <directory of recorded json>
  python json_expect.py back   <out-dir> <directory of recorded json>

For each case the C++ test writes <out>/<name>.json in the same tagged form the
generator used. Comparing whole trees is the point: checking a handful of keys
by hand let a broken reader pass (7th review).
"""
import json
import math
from pathlib import Path
import struct
import sys


def bits(value):
    return format(struct.unpack('<Q', struct.pack('<d', value))[0], '016x')


def marked(value):
    """The same tagging the JSON cases use for a value a reader must refuse."""
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, int):
        return dict(int=str(value))
    if isinstance(value, float):
        return dict(float=bits(value))
    if isinstance(value, str):
        # The bytes travel with the text: a reader that keeps a \\u0000 escape as
        # six ASCII characters and a writer that fails to escape the backslash
        # cancel out if only the re-encoded text is compared (8th review).
        return dict(str=value, bytes=value.encode('utf-8', 'surrogatepass').hex())
    if isinstance(value, list):
        return [marked(v) for v in value]
    # The key is tagged like any other string, bytes included: a reader that kept
    # "\\u0061" as six characters for a KEY passed while only values were checked
    # (9th review), and at("a") would then fail on real input.
    return dict(object=[[marked(k), marked(v)] for k, v in value.items()])


def cleaned(value):
    """What adaptive_refine.clean() does before saving: a non-finite float is a string.

    The writer follows that rule, so a round trip turns NaN and Infinity into
    "nan"/"inf"/"-inf" by design. Without this the comparison reports thousands of
    differences that are the convention working (measured on the fuzz corpus).
    """
    if isinstance(value, float) and not math.isfinite(value):
        return {float('inf'): 'inf', float('-inf'): '-inf'}.get(value, 'nan')
    if isinstance(value, list):
        return [cleaned(v) for v in value]
    if isinstance(value, dict):
        return {k: cleaned(v) for k, v in value.items()}
    return value


def safe(text):
    """Printable on a UTF-8 stream: a lone surrogate cannot be encoded."""
    return text.encode('utf-8', 'backslashreplace').decode()


def nan(text):
    """True for any binary64 NaN pattern: NaN equals NaN (design 3.3)."""
    value = int(text, 16)
    return (value >> 52) & 0x7ff == 0x7ff and value & ((1 << 52)-1)


def differences(expected, found, path='/'):
    # Type first, always: Python's == says True equals 1, so a reader that turned
    # true into 1 compared equal to its expectation (8th review).
    if type(expected) is not type(found):
        return [f'{path}: {type(expected).__name__} {expected!r} vs {type(found).__name__} {found!r}']
    if isinstance(expected, dict) and isinstance(found, dict):
        if set(expected) != set(found):
            return [f'{path}: tags {sorted(expected)} vs {sorted(found)}']
        out = []
        for key in expected:
            if key == 'object':
                left, right = expected[key], found[key]
                if len(left) != len(right):
                    out.append(f'{path}: {len(left)} keys vs {len(right)}')
                    continue
                for i, ((lk, lv), (rk, rv)) in enumerate(zip(left, right)):
                    # Keys are compared in place, so their order, text and bytes all count.
                    name = lk.get('str') if isinstance(lk, dict) else lk
                    out += differences(lk, rk, f'{path}{name}=key/')
                    out += differences(lv, rv, f'{path}{name}/')
            elif key == 'float':
                if not (nan(expected[key]) and nan(found[key])) and expected[key] != found[key]:
                    out.append(f'{path}: float {expected[key]} vs {found[key]}')
            else:
                out += differences(expected[key], found[key], path)
        return out
    if isinstance(expected, list) and isinstance(found, list):
        if len(expected) != len(found):
            return [f'{path}: {len(expected)} items vs {len(found)}']
        out = []
        for i, (a, b) in enumerate(zip(expected, found)):
            out += differences(a, b, f'{path}{i}/')
        return out
    return [] if expected == found else [f'{path}: {expected!r} vs {found!r}']


def flat(value, out=None):
    """The same content as a walk: object/array boundaries, keys and tagged leaves.

    The C++ side writes this form for documents whose tagged tree would be deeper
    than the supported nesting. Comparing the walk keeps those files verified
    instead of set aside (11th review).
    """
    out = [] if out is None else out
    if isinstance(value, dict):
        out.append({'o': '{'})
        for key, item in value.items():
            out.append({'key': marked(key)})
            flat(item, out)
        out.append({'o': '}'})
    elif isinstance(value, list):
        out.append({'a': '['})
        for item in value:
            flat(item, out)
        out.append({'a': ']'})
    else:
        out.append(marked(value))
    return out


INT64 = range(-2**63, 2**63)


def read_with_limit(path):
    """CPython reading a document under this port's limit, as an expectation.

    CPython's integers are unbounded and this port's are int64 (plan stage 7), so
    a document CPython reads and this reader refuses is not a defect - it is the
    declared limit. json.loads calls parse_int on each integer token in document
    order, which is where this reader stops, so the refusal lands in the same
    place with the same message.
    """
    def within(digits):
        value = int(digits)
        if value not in INT64:
            raise ValueError(f'integer {digits} does not fit in 64 bits')
        return value
    try:
        value = json.loads(path.read_bytes().decode('utf-8'), parse_int=within)
        return dict(ok=True, value=marked(value), raw=value)
    except ValueError as error:
        return dict(ok=False, error=[type(error).__name__, str(error)],
                    pos=getattr(error, 'pos', None))


def compare(name, expected, found):
    if expected['ok'] and not found['ok'] and found['error'][0] == 'RecursionError':
        # The declared nesting limit, which is narrower than CPython's on purpose
        # (design 3.2). Anything else that only one side accepts is a mismatch.
        return []
    if expected['ok'] != found['ok']:
        return [f"{name}: ok={found['ok']}, expected {expected['ok']}"]
    if expected['ok']:
        if 'flat' in found:
            if 'raw' not in expected:
                return [f'{name}: written as a walk, but there is nothing to walk against']
            return [f'{name} {line}' for line in differences(flat(expected['raw']), found['flat'])]
        return [f'{name} {line}' for line in differences(expected['value'], found['value'])]
    problems = []
    if expected['error'][0] != found['error'][0]:
        problems.append(f"{name}: {found['error'][0]}, expected {expected['error'][0]}")
    if expected['error'][1] != found['error'][1]:
        problems.append(f"{name}: message {found['error'][1]!r}, expected {expected['error'][1]!r}")
    if expected.get('pos') != found.get('pos'):
        problems.append(f"{name}: pos {found.get('pos')}, expected {expected.get('pos')}")
    return problems


def main(argv):
    # The tool walks structures as deep as the port supports (500), and marked()
    # and flat() add frames of their own; Python's default 1,000 is not enough.
    # CPython's json refuses about 9,999 levels, so nothing deeper can arrive.
    sys.setrecursionlimit(50000)
    mode, out, source = argv[0], Path(argv[1]), Path(argv[2])
    problems, checked = [], 0
    # The byte comparison is only asked for where Python's own writer made the file.
    want_bytes = len(argv) > 3 and argv[3] == 'bytes'
    bytes_same = []
    if mode == 'back':
        # What our writer produced must read back as the very same value: same types,
        # same key order, same floats to the bit.
        for path in sorted(out.rglob('*.json')):
            name = path.relative_to(out).as_posix()
            original = source/name
            if not original.is_file():
                problems.append(f'{name}: written back from nowhere')
                continue
            checked += 1
            problems += [f'{name} {line}' for line in
                         differences(marked(cleaned(json.loads(original.read_text()))),
                                     marked(json.loads(path.read_text())))]
            # Stronger than equal values, where it is asked for: byte for byte. Only
            # files written with python3's json.dumps(indent=2) qualify - the
            # response.* files are the CCEX binary's, and a fuzz corpus is written by
            # the generator, both with other conventions.
            if not want_bytes or path.name.startswith('response.'):
                pass
            elif path.read_bytes() != original.read_bytes():
                problems.append(f'{name}: same value, different bytes')
            else:
                bytes_same.append(True)
        for problem in problems:
            print('MISMATCH', safe(problem))
        print(f'{checked} files compared, {len(bytes_same)} of them byte for byte')
        return 1 if problems or not checked or (want_bytes and not bytes_same) else 0
    if mode == 'cases':
        wanted = []
        for p in sorted(source.glob('*.read.json')):
            name = p.name[:-len('.read.json')]
            expected = json.loads(p.read_text())
            if expected['ok']:
                # The stored expectation is CPython's, whose integers are wider
                # than this port's; where that is what separates them, the port
                # has to refuse rather than approximate.
                again = read_with_limit(source/f'{name}.json')
                if not again['ok']:
                    expected = again
            wanted.append((name, expected))
    else:
        wanted = [(path.relative_to(source).as_posix()[:-len('.json')], read_with_limit(path))
                  for path in sorted(source.rglob('*.json'))]
    for name, expected in wanted:
        found_path = out/f'{name}.json'
        if not found_path.is_file():
            problems.append(f'{name}: the C++ reader wrote nothing')
            continue
        checked += 1
        problems += compare(name, expected, json.loads(found_path.read_text()))
    names = {name for name, _ in wanted}
    extra = sorted(p.relative_to(out).as_posix()[:-len('.json')] for p in out.rglob('*.json'))
    extra = [name for name in extra if name not in names]
    if extra:
        # A leftover output from an earlier run would otherwise be compared and pass
        # while the new reader wrote nothing (8th review).
        problems.append(f'unexpected files under {out}: {extra[:5]}')
    for problem in problems:
        print('MISMATCH', safe(problem))
    print(f'{checked} files compared')
    return 1 if problems or not checked else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
