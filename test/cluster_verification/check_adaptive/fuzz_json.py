#!/usr/bin/env python3
"""Random JSON, read by CPython, for the C++ reader to be compared against (plan 0B-1).

  python fuzz_json.py OUT_DIR COUNT SEED

Writes COUNT files OUT_DIR/f<NNNNN>.json. Roughly half are malformed: the value
of a differential test is in the inputs nobody thought to list, and every defect
review found in this reader was of that shape - the fixed case list passed while
the implementation was wrong.

The comparison is the usual pair:
  ./bin/test_json_read --all OUT_DIR
  python json_expect.py tree <out>/all OUT_DIR
"""
import json
import random
import sys
from pathlib import Path

ESCAPES = ['\\"', '\\\\', '\\/', '\\b', '\\f', '\\n', '\\r', '\\t',
           '\\u0000', '\\u0041', '\\u00e9', '\\uffff', '\\ud800', '\\udbff\\udfff',
           '\\ud800\\ud801\\udc00', '\\udc00', '\\ud83d\\ude00']
NUMBERS = ['0', '-0', '1', '-1', '0.0', '-0.0', '1e0', '1E+2', '1e-400', '1e400', '2.0',
           '9007199254740993', '-9223372036854775809', '18446744073709551616',
           '0.1', '123.456e-7', '1' + '0'*30, '-' + '9'*40]
BROKEN = ['01', '1e', '+1', '.5', '5.', '--1', '0x10', 'nan', 'None', "'x'", '1.2.3',
          '"\\q"', '"\\u12"', '"unterminated', 'tru', '{,}', '[,]', '{"a" 1}', '{"a": }',
          '{"a": 1,}', '[1,]', '{1: 2}', '', '   ', '}', ']', '{"a": 1} 2', '[1] [2]']
TEXTS = ['', 'x', 'key', 'a b', '0', 'true', '\\u0073eed', 'é', '😀', ' spaced ']


def text(rng):
    parts = [rng.choice(TEXTS)]
    for _ in range(rng.randint(0, 3)):
        parts.append(rng.choice(ESCAPES) if rng.random() < .6 else rng.choice(TEXTS))
    return '"'+''.join(parts)+'"'


def value(rng, depth):
    pick = rng.random()
    if depth > 4 or pick < .25:
        return rng.choice([rng.choice(NUMBERS), text(rng), 'true', 'false', 'null',
                           'NaN', 'Infinity', '-Infinity'])
    space = rng.choice(['', ' ', '\n', '\t', '  '])
    if pick < .6:
        items = [value(rng, depth+1) for _ in range(rng.randint(0, 4))]
        return '['+space+(','+space).join(items)+space+']'
    pairs = [f'{text(rng)}{space}:{space}{value(rng, depth+1)}' for _ in range(rng.randint(0, 4))]
    return '{'+space+(','+space).join(pairs)+space+'}'


def damaged(rng, body):
    """One mutation of a well-formed document, or an outright broken token."""
    how = rng.randint(0, 3)
    if how == 0 or not body:
        return rng.choice(BROKEN)
    if how == 1:
        return body[:rng.randint(0, max(0, len(body)-1))]          # truncate
    if how == 2:
        cut = rng.randrange(len(body))
        return body[:cut]+rng.choice(',:{}[]"\\ e+-01')+body[cut+1:]   # replace one character
    cut = rng.randrange(len(body)+1)
    return body[:cut]+rng.choice([',', ':', '"', '\\', ']', '}'])+body[cut:]   # insert one


def main(argv):
    out, count, seed = Path(argv[0]), int(argv[1]), int(argv[2])
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    accepted = 0
    for i in range(count):
        body = value(rng, 0)
        if rng.random() < .5:
            body = damaged(rng, body)
        path = out/f'f{i:05d}.json'
        path.write_text(body)
        try:
            json.loads(body)
            accepted += 1
        except ValueError:
            pass
    print(f'{count} files, CPython accepts {accepted}, refuses {count-accepted}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
