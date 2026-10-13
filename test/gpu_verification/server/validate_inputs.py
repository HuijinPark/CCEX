#!/usr/bin/env python3
"""Guard the server benchmark's authorized order-two, unmodified bath scope."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
INPUTS = ROOT / 'test/code_verification/check_input_forms/inputs'


def main():
    checked = []
    for case in sorted(INPUTS.glob('Diamond_*')):
        p1 = '_P1_' in case.name
        expected = dict(order=2, rbath=200 if p1 else 30, rdip=100 if p1 else 8,
                        nstep=40 if p1 else 30, pulse={'count': 1}, seed=1)
        for state in (0, 1):
            path = case / 'inputfiles' / f'ccein_nstate{state}.json'
            config = json.loads(path.read_text())
            for key, value in dict(expected, nstate=state).items():
                if config.get(key) != value:
                    raise ValueError(f'{path}: expected {key}={value}, got {config.get(key)}')
            method = 'gCCE' if case.name.endswith('_gcce') else 'cCCE'
            if config.get('method') != method:
                raise ValueError(f'{path}: unexpected method')
            checked.append(dict(case=case.name, nstate=state, **expected))
    if len(checked) != 8:
        raise ValueError(f'Expected eight cases, found {len(checked)}')
    print(json.dumps(checked, indent=2))


if __name__ == '__main__':
    main()

