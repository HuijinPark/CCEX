"""Bath-state files for the NV-13C example: for every bath file config_<k> in DATA (from BathGeneratorSpinCluster), NSTATE
files cfg<k>_state<s> (+-1/2 per spin, in the order of the bath file lines) with default_rng(20261007 + k).
    make_states.py DATA [NSTATE]"""
import re, sys
from pathlib import Path
import numpy as np

data = Path(sys.argv[1]); nstate = int(sys.argv[2]) if len(sys.argv) > 2 else 30
files = sorted((p for p in data.glob("config_*") if re.fullmatch(r"config_\d+", p.name)), key=lambda p: int(p.name[7:]))
for f in files:
    k = int(f.name[7:])
    n = int(round(float(f.read_text().split()[0]))) - 1
    rng = np.random.default_rng(20261007 + k)
    for s in range(1, nstate + 1):
        (data / f"cfg{k}_state{s}").write_text(f"{n}\n" + "".join(f"{v:.5f}\n" for v in rng.choice([0.5, -0.5], n)))
print(f"{len(files)} configurations x {nstate} states -> {data}")
