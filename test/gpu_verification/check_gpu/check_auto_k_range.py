#!/usr/bin/env python3
"""backend auto with configuration.k_range on two ranks, when only one rank opens the card (A73).

From the second configuration on, the rank that had opened the card for an eligible cluster skipped the collective node
split (MPI_Comm_split_type) that the other rank entered, and the run hung (review 20261015 feature-dev 1). The bath:
one triple of 13C (dimension 8, the only cluster at gpu_options.min_dim 8, so one rank gets it), a pair and three
singles; cCCE3 Hahn, k_range [1, 2] over two copies of that bath. Pass: exit 0 within the time limit, exactly one rank
opened the card, both configurations written and the same bytes (the same bath).

Usage: python3 check_auto_k_range.py --cuda bin/main_gpu.out --output DIR [--mpirun MPIRUN]
"""
import argparse
import json
import os
import shutil
import signal
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cuda', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mpirun', default=shutil.which('mpirun') or 'mpirun')
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    rows = [(5.0, 0, 0), (7.0, 0, 0), (6.0, 1.7, 0),     # the triple, about 2 A apart
            (0, 20, 0), (0, 21.5, 0),                     # a pair
            (-20, 0, 0), (0, 0, 20), (0, -20, 0)]         # singles
    for k in (1, 2):
        lines = [f"{len(rows)}.00000\t0.00000\t0.00000\t0.00000"] + [f"{x:.5f}\t{y:.5f}\t{z:.5f}\t13C" for x, y, z in rows]
        (out / f"bath_{k}").write_text("\n".join(lines) + "\n")
    config = {
        "method": "cCCE", "quantity": "coherence",
        "Qubit": {"qubit": [{"name": "NV", "spin": 1.0, "xyz": [0.0, 0.0, 0.0]}], "alpha_state": {"ms": [1]}, "beta_state": {"ms": [0]}},
        "configuration": {"bathfile": ["bath_{k}"], "k_range": [1, 2], "seed": 1},
        "order": 3, "bfield": 500, "rbath": 30, "rdip": 3, "deltat": 0.005, "nstep": 50, "pulse": {"count": 1},
        "execution": {"backend": "auto", "gpu_options": {"min_dim": 8, "packed_max_dim": 0, "warp_max_dim": 0}},
        "output": {"file": "auto_k"}}
    (out / "input.json").write_text(json.dumps(config, indent=1))
    env = dict(os.environ, I_MPI_FABRICS='shm')
    with (out / "run.log").open("w") as log:
        # a session of its own: on a hang the ranks go with mpirun (killing mpirun alone left them running)
        p = subprocess.Popen([args.mpirun, '-n', '2', str(args.cuda.resolve()), '-f', 'input.json'], cwd=out, env=env,
                             stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True)
        try:
            p.wait(timeout=180)
        except subprocess.TimeoutExpired:
            os.killpg(p.pid, signal.SIGKILL); p.wait()
            raise SystemExit(f"FAIL auto-k-range: no end in 180 s (the hang of the node split; see {out}/run.log)")
    text = (out / "run.log").read_text()
    assert p.returncode == 0, ("exit", p.returncode, text[-1500:])
    opened = [l for l in text.splitlines() if l.startswith("[backend] rank=")]
    assert len(opened) == 1, ("the case needs exactly one rank to open the card", opened)
    for name in ("wiDiv", "noDiv"):
        a, b = out / f"auto_k_k1_{name}", out / f"auto_k_k2_{name}"
        assert a.is_file() and b.is_file() and a.read_bytes() == b.read_bytes(), name
    print("PASS auto-k-range (two ranks, one opened the card, both configurations, the same bytes)", flush=True)


if __name__ == '__main__':
    main()
