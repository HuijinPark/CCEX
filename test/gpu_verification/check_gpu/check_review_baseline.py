#!/usr/bin/env python3
"""Replay completed CPU integration inputs with a frozen pre-review binary."""
import argparse
import json
from pathlib import Path
from benchmark import compare
from check_pcce import run_case


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--results', type=Path, required=True)
    args = p.parse_args()
    records = []
    for suite in ('acceleration', 'pulses', 'pcce'):
        for source in sorted((args.results/suite).glob('*/cpu/input.json')):
            config = json.loads(source.read_text())
            directory = args.results/'baseline'/suite/source.parent.parent.name
            execution, log = run_case(config, directory.resolve(), args.baseline.resolve(), ranks=2)
            assert execution['returncode'] == 0, log[-2000:]
            error = compare(source.parent, directory, 1e-10, 1e-10)
            records.append(dict(case=f'{suite}/{source.parent.parent.name}', **error))
            print('PASS baseline', records[-1], flush=True)
    assert len(records) == 28, len(records)
    (args.results/'baseline/summary.json').write_text(json.dumps(dict(status='passed', cases=records), indent=2)+'\n')


if __name__ == '__main__':
    main()
