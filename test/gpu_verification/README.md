# test/gpu_verification/

The CUDA backend (`src/calculation/gpu/`, `backend` in the input). Three folders.

| folder | what | needs |
|---|---|---|
| `check_gpu/` | the checks: 14 `test_*.cpp` unit tests (which `make` target builds and runs which: `check_gpu/README.md`), nine `check_*.py` that run CPU and CUDA builds on the same inputs and compare (pulses, gCCE4 accuracy and batching, pCCE, seeds and sampled states, options, host leaks); `examples/` and `data/` are the small inputs they use; the seven `validation-*.json` are the recorded results (internal notes) cites | a CUDA build (`bin/main_gpu.out`, idm4-gpu) for most; `check_review_controls.py --cpu` and `check_seed_states.py --cpu` run on the CPU build |
| `gpu_benchmarks/` | timing, not correctness: CPU against CUDA, rank scaling, the PL6 pCCE case | idm4-gpu |
| `server/` | the three standing Slurm jobs for idm4-gpu (build, CCE2 measurement, high orders) with the environment they use, last run on 2026-10-02 (jobs 640, 641, 644), and two measurement jobs (`min_dim`, `gpu_sharing`); `README.md` there. `make gpu-test` passed in idm4-gpu jobs 691 (83/0) and 692 (84/0) on 2026-10-12 | idm4-gpu |
