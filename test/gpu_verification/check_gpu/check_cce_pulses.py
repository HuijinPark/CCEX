#!/usr/bin/env python3
"""Compare CPU/GPU cCCE Ramsey, odd/even CPMG, and manual sampled intervals."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

from benchmark import INPUTS, compare, run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpu',type=Path,required=True)
    parser.add_argument('--cuda',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parsed=parser.parse_args()
    root=parsed.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    args=SimpleNamespace(nearest_sites=3,steps=9,batch=4,min_dim=None,evolution=None,
                         mpirun=None,ranks=1,timeout=120,savemode='coherence',pulses=None)
    source_root=INPUTS/'Diamond_NV_P1_10ppm_cce'
    records=[]

    def check(name,source):
        directory=root/name
        cpu=run(parsed.cpu.resolve(),source,directory/'cpu','cpu',args)
        gpu=run(parsed.cuda.resolve(),source,directory/'cuda','cuda',args)
        error=compare(directory/'cpu',directory/'cuda',1e-8,2e-6)
        records.append(dict(case=name,reference=cpu,target=gpu,**error))
        (root/'summary.json').write_text(json.dumps(records,indent=2)+'\n')
        print(f'PASS {name}: max_error={error["max_abs_error"]:.3g}',flush=True)

    for state in [0,1]:
        source=source_root/'inputfiles'/f'ccein_nstate{state}.json'
        for pulses in [0,2,3,4]:
            args.pulses=pulses
            check(f'nstate{state}_pulses{pulses}',source)

    # Keep the fixture's absolute file paths while changing only the sequence.
    config=json.loads((root/'nstate1_pulses3/cpu/input.json').read_text())
    config.update(pulsename='manual',sequence=[['1/8','X',180],['2/8','Y',180],
                                             ['3/8','X',180],['2/8','I',0]])
    source=root/'manual.json'
    source.write_text(json.dumps(config))
    args.pulses=None
    check('manual_single',source)


if __name__=='__main__':
    main()
