"""Bath-state files for the NV-13C example: for every bath file config_<k> in DATA (make_data.sh), NSTATE files
cfg<k>_state<s> (m = +-1/2 per 13C, in the order of the bath file lines) with the rule of example/tools/StateFiles
(numpy default_rng([seed, k, s]), seed 1 as in input.jsonc; the rule the runner uses for the states it makes, A94 - until
2026-10-16 this script drew them with default_rng(20261007 + k), so the files and the example's results changed then).
    make_states.py DATA [NSTATE]"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "StateFiles"))
from make_state_files import write_bath_states  # noqa: E402

SEED = 1                                  # configuration.seed of input.jsonc
data = Path(sys.argv[1]); nstate = int(sys.argv[2]) if len(sys.argv) > 2 else 30
files = sorted((p for p in data.glob("config_*") if re.fullmatch(r"config_\d+", p.name)), key=lambda p: int(p.name[7:]))
for f in files:
    k = int(f.name[7:])
    for s in range(1, nstate + 1):
        write_bath_states(data / f"cfg{k}_state{s}", f, {"13C": 0.5}, SEED, k, s)
print(f"{len(files)} configurations x {nstate} states -> {data}")
