#!/usr/bin/env python3
"""Judge the config reader's verdicts with CPython's parse_int hook (plan 0B-1).

  python config_verdicts.py CORPUS_DIR VERDICTS_TSV

`test_json_config --verdicts` prints one line per file: "ok" or the exception.
CPython's json has a parse_int hook that fires on every integer TOKEN, duplicates
included, which is exactly what the config check claims to do - so the two can be
compared without writing a second parser. The path in the message differs (the
hook does not know it), so the comparison is the verdict, and the integer named.
"""
import json
import re
import sys
from pathlib import Path

LIMIT = 2**63


class Rejected(Exception):
    def __init__(self, value):
        self.value = value


def by_python(path):
    try:
        text = path.read_bytes().decode('utf-8')
    except UnicodeDecodeError:
        return 'error'

    def check(text):
        value = int(text)
        if not -LIMIT <= value < LIMIT:
            raise Rejected(value)
        return value

    try:
        json.loads(text, parse_int=check)
        return 'ok'
    except Rejected as rejected:
        return f'reject {rejected.value}'
    except (ValueError, RecursionError):
        return 'error'


def by_cpp(verdict):
    if verdict == 'ok':
        return 'ok'
    found = re.match(r'ValueError: integer (-?\d+) at ', verdict)
    return f'reject {found.group(1)}' if found else 'error'


def main(argv):
    corpus, verdicts = Path(argv[0]), Path(argv[1])
    problems, total = [], 0
    for line in verdicts.read_bytes().decode('utf-8', 'replace').splitlines():
        if '\t' not in line:
            continue
        name, verdict = line.split('\t', 1)
        total += 1
        ours, theirs = by_cpp(verdict), by_python(corpus/name)
        if ours != theirs:
            problems.append(f'{name}: {ours}, expected {theirs} ({verdict[:60]})')
    for problem in problems[:20]:
        print('MISMATCH', problem)
    print(f'{total} verdicts compared')
    return 1 if problems or not total else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
