"""The input reading of the automatic convergence (example/tools/AutoConvergence): a CCEX input with its comments (//, #
and !, outside strings, as CCEX reads them), its keys in the names the tools look up (case not read, as in CCEX), the
A67 units, and the Settings `python -m autoconv analyze` and `judge` take from the input's auto_convergence block - the
runner's block and keys (A94, G19; the driver configuration of its own, {"base_input", "auto_convergence": {"tolerance",
"steps", ...}, "bath_source", "state_source"}, was removed with its key system).
"""
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

# The error budget of the analysis (section 0.4 of the plan): truncation and numerics 1/3 of the total, sampling 2/3; inside
# the first third rdip : rbath = 50 : 50. The analysis's own settings, reported beside its estimates (not input keys).
BUDGET_TRUNCATION, BUDGET_SAMPLING = 1 / 3, 2 / 3
SPLIT = {"rdip": 0.5, "rbath": 0.5}


class ConfigError(ValueError):
    """An input (its auto_convergence block) the analysis refuses, with the reason."""


def strip_comments(text):
    """Remove '#', '!' and '//' comments (all accepted by the CCEX input reader) outside JSON strings."""
    out, in_str, esc, skip = [], False, False, False
    for i, ch in enumerate(text):
        if skip:
            if ch == "\n":
                skip = False
                out.append(ch)
            continue
        if in_str:
            out.append(ch)
            if esc:
                esc = False
            elif ch == "\\":
                esc = True
            elif ch == '"':
                in_str = False
        elif ch == '"':
            in_str = True
            out.append(ch)
        elif ch in "#!" or text.startswith("//", i):
            skip = True
        else:
            out.append(ch)
    return "".join(out)


def read_json(path):
    """JSON with '#', '!' or '//' comments (CCEX input style)."""
    return json.loads(strip_comments(Path(path).read_text()))


# CCEX reads input keys without regard to case (docs/user/InputFile.md) and the runner looks them up by their written
# names, so a CCEX input's keys are put in those names when it is read (review 20261013-180049 R4): lower case, but for
# the few names written with capitals. Values (paths, names, strings) stay as written.
MIXED_CASE_KEYS = {"ccce": "cCCE", "gcce": "gCCE", "mecce": "meCCE", "delta_a": "delta_A", "d": "D", "e": "E"}
TOP_MIXED_CASE_KEYS = {"qubit": "Qubit", "defect": "Defect"}      # Qubit.qubit stays lower case


class KeyCaseError(ValueError):
    """Two keys of one block that differ only in case."""


def _single_keys(pairs):
    seen = {}
    for k, _ in pairs:
        if k.lower() in seen:
            raise KeyCaseError(f'input keys "{seen[k.lower()]}" and "{k}" are one key to CCEX (case is not read); keep one')
        seen[k.lower()] = k
    return dict(pairs)


def canonical_keys(value, top=True):
    if isinstance(value, list):
        return [canonical_keys(v, False) for v in value]
    if not isinstance(value, dict):
        return value
    return {(TOP_MIXED_CASE_KEYS.get(k.lower()) if top else None) or MIXED_CASE_KEYS.get(k.lower(), k.lower()):
            canonical_keys(v, False) for k, v in value.items()}


class UnitError(ValueError):
    """A key written as {"value", "unit"} that CCEX would refuse (A67)."""


# A67: the keys CCEX reads in a unit of their own, written as a plain number (that unit) or {"value": number or list,
# "unit": ...}. The same table and the same exact operation as the CCEX reader (src/input/option.cpp, applyInputUnits): a
# multiplication by an exact power of ten towards a smaller unit, a division by one towards a larger, so the inputs the
# runner writes hold the number CCEX itself would have read.
UNITS = {"ms": {"s": (1e3, 1), "ms": (1, 1), "us": (1, 1e3), "ns": (1, 1e6)},
         "ns": {"s": (1e9, 1), "ms": (1e6, 1), "us": (1e3, 1), "ns": (1, 1)},
         "G": {"G": (1, 1), "mT": (10, 1), "T": (1e4, 1)},
         "A": {"A": (1, 1), "nm": (10, 1)},
         "MHz": {"Hz": (1, 1e6), "kHz": (1, 1e3), "MHz": (1, 1), "GHz": (1e3, 1)},
         "kHz": {"Hz": (1, 1e3), "kHz": (1, 1), "MHz": (1e3, 1), "GHz": (1e6, 1)},
         # A94 (G20): cluster.bundle.merge_cutoff in rad/ms
         "rad/ms": {"rad/s": (1, 1e3), "rad/ms": (1, 1), "rad/us": (1e3, 1)},
         # A94 (G20): Qubit.qubit[].gyro in rad kHz/G with the table of the Defect gyros (option.cpp,
         # readDefectGyroScaleToRadkHzPerG: value x the factor, not a power of ten)
         "radkHz/G": {"radkHz/G": (1, 1), "rad/ms/G": (1, 1), "kHz/G": (2.0 * math.pi, 1), "MHz/G": (2.0 * math.pi * 1.0e3, 1),
                      "MHz/T": (2.0 * math.pi * 1.0e-1, 1), "Hz/T": (2.0 * math.pi * 1.0e-7, 1), "radMHz/T": (1.0e-1, 1),
                      "rad/us/T": (1.0e-1, 1), "rad/s/T": (1.0e-7, 1)}}


def plain_value(where, value, own):
    """value in the unit `own` (a key of UNITS): a plain value as it is, {"value", "unit"} converted; UnitError otherwise."""
    if not isinstance(value, dict):
        return value
    units = UNITS[own]
    keys = {k.lower(): v for k, v in value.items()}

    def numbers(v):
        return (isinstance(v, (int, float)) and not isinstance(v, bool)) or (isinstance(v, list) and v and all(map(numbers, v)))
    if set(keys) != {"value", "unit"} or not isinstance(keys["unit"], str) or not numbers(keys["value"]):
        raise UnitError(f'{where} must be a number (in {own}) or {{"value": number or list of numbers, "unit": "..."}} with '
                        f'the unit one of {", ".join(units)}')
    if keys["unit"] not in units:
        raise UnitError(f'{where}.unit "{keys["unit"]}" is not known; the unit is one of {", ".join(units)} (a plain number is in {own})')
    times, over = units[keys["unit"]]

    def scale(v):
        if isinstance(v, list):
            return [scale(x) for x in v]
        return v * times if times != 1 else v / over if over != 1 else v
    return scale(keys["value"])


def _get(d, *path):
    for k in path:
        d = d.get(k) if isinstance(d, dict) else None
    return d


def plain_units(d):
    """The input d (read_input's names) with every A67 key in its own unit, in place."""
    def at(block, key, where, own):
        if isinstance(block, dict) and key in block:
            block[key] = plain_value(where, block[key], own)
    for key, own in (("deltat", "ms"), ("bfield", "G"), ("rbath", "A"), ("rdip", "A")):
        at(d, key, key, own)
    at(d.get("configuration"), "bath_shift", "configuration.bath_shift", "A")
    qubits = d.get("Qubit", {}).get("qubit") if isinstance(d.get("Qubit"), dict) else None
    for i, q in enumerate(qubits if isinstance(qubits, list) else []):
        at(q, "xyz", f"Qubit.qubit[{i}].xyz", "A")
        at(q, "detuning", f"Qubit.qubit[{i}].detuning", "kHz")
        at(q, "gyro", f"Qubit.qubit[{i}].gyro", "radkHz/G")                      # A94 (G20)
    at(_get(d, "cluster", "bundle"), "merge_cutoff", "cluster.bundle.merge_cutoff", "rad/ms")
    at(_get(d, "hamiltonian", "common", "quadrupole"), "cell", "hamiltonian.common.quadrupole.cell", "A")
    pulse = d.get("pulse")
    if isinstance(pulse, dict):
        at(pulse, "duration", "pulse.duration", "ns")
        for key in ("energy_shift", "bandwidth"):
            at(pulse.get("bath"), key, f"pulse.bath.{key}", "MHz")
    ham = d.get("hamiltonian")
    for name, section in (ham.items() if isinstance(ham, dict) else ()):
        jumps = section.get("jump_operators") if isinstance(section, dict) else None
        for i, j in enumerate(jumps if isinstance(jumps, list) else []):
            at(j, "rate", f"hamiltonian.{name}.jump_operators[{i}].rate", "MHz")
    return d


def read_input(path):
    """A CCEX input (read_json) with its keys in the names the tools look up; keys that differ only in case are refused."""
    return canonical_keys(json.loads(strip_comments(Path(path).read_text()), object_pairs_hook=_single_keys))


def _fraction(where, value, lo_open=0.0, hi_open=1.0):
    if not (isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value) and lo_open < value < hi_open):
        raise ConfigError(f"auto_convergence.{where} = {value!r} must be a number between {lo_open:g} and {hi_open:g}")
    return float(value)


@dataclass
class Settings:
    """What `python -m autoconv analyze` and `judge` use, from the auto_convergence block of a CCEX input (A94, G19: the
    runner's block; the driver configuration of its own, with keys the runner refuses, was removed). steps: the parts the
    analysis judges - the time window unless time_window.value is given, the sampling unless sampling.configurations is
    (rdip and rbath are the runner's, never judged here). The error budget (tolerance and its split), the error norm, the
    window rule and the truncation tolerances of judge are the analysis's own fixed settings, not input keys."""
    steps: dict
    decay_threshold: float = 0.01          # time_window.decay_level
    window_margin: float = 1.5             # time_window.production_factor
    t2_fit_points: str = "peaks"           # t2_fit.points
    t2_se_target: float = 0.05             # sampling.threshold (the SE(A)/A target, as the runner's)
    tolerance: float = 0.01
    rdip_tolerance: float = 0.01 * BUDGET_TRUNCATION * SPLIT["rdip"]
    rbath_tolerance: float = 0.01 * BUDGET_TRUNCATION * SPLIT["rbath"]
    sampling_tolerance: float = 0.01 * BUDGET_SAMPLING
    error_norm: str = "area"
    window_method: str = "fit"             # T_d where the T2 fit line reaches decay_threshold ("envelope": the envelope itself)
    window_median_points: int = 9
    truncation_max_tol: float = 0.01
    truncation_area_tol: float = 0.005
    base_input: dict = field(default_factory=dict)
    root: Path = Path(".")


def parse_auto_convergence(block=None, base_input=None, root="."):
    """Settings from an auto_convergence block (the runner's keys; old or unknown keys refused as the runner refuses them,
    ConfigError). base_input: the rest of the CCEX input (for the expected run shape and the bath files)."""
    from .runner import Refused, Runner, WHERE         # the runner's table of the block's keys (one key system)
    block = {} if block is None else block
    if not isinstance(block, dict):
        raise ConfigError("auto_convergence must be an object")
    try:
        flat = Runner._flatten({k: v for k, v in block.items() if k != "enabled"})
    except Refused as e:
        raise ConfigError(str(e))
    s = Settings(steps=dict(time_window="set_time_window" not in flat, rdip=False, rbath=False,
                            sampling="number_config" not in flat))
    if "threshold_time_window" in flat:
        s.decay_threshold = _fraction(WHERE["threshold_time_window"], flat["threshold_time_window"])
    if "production_window_factor" in flat:
        v = flat["production_window_factor"]
        if not (isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) and v >= 1):
            raise ConfigError(f"auto_convergence.time_window.production_factor = {v!r} must be at least 1")
        s.window_margin = float(v)
    if "fit_points" in flat:
        if flat["fit_points"] not in ("peaks", "all"):
            raise ConfigError(f"auto_convergence.t2_fit.points must be peaks or all, not {flat['fit_points']!r}")
        s.t2_fit_points = flat["fit_points"]
    if "threshold_config" in flat:
        s.t2_se_target = _fraction(WHERE["threshold_config"], flat["threshold_config"])
    s.base_input, s.root = dict(base_input or {}), Path(root)
    return s


def load_config(path):
    """Settings of a CCEX input with an auto_convergence block (the input `main.out -f` runs): its block, and the rest as
    base_input. Comments //, # and ! as CCEX reads them."""
    path = Path(path)
    try:
        cin = read_input(path)
        plain_units(cin)
    except (KeyCaseError, UnitError) as e:
        raise ConfigError(str(e))
    if "auto_convergence" not in cin:
        raise ConfigError(f"{path}: the input has no auto_convergence block")
    block = cin.pop("auto_convergence")
    return parse_auto_convergence(block, cin, path.parent)


def source_capacity(s):
    """Bath files of the input's configuration.bathfile pattern ({k}) that exist: the k of configuration.k_range with a
    file, or (without k_range) every k found; None without a pattern."""
    from . import autostart
    conf = s.base_input.get("configuration", {}) if isinstance(s.base_input, dict) else {}
    names = conf.get("bathfile") if isinstance(conf, dict) else None
    if not (isinstance(names, list) and len(names) == 1 and isinstance(names[0], str) and "{k}" in names[0]):
        return None
    pattern = names[0] if conf.get("path") is None or Path(names[0]).is_absolute() else str(Path(conf["path"]) / names[0])
    kr = conf.get("k_range")
    if isinstance(kr, list) and len(kr) == 2:
        p = lambda k: Path(pattern.replace("{k}", str(k)))
        return sum(1 for k in range(kr[0], kr[1] + 1) if (p(k) if p(k).is_absolute() else Path(s.root) / p(k)).is_file())
    try:
        return len(autostart.discover(pattern, s.root))
    except ValueError:
        return None


def read_manifest(path):
    """Pilot/production manifest: {"configurations": [{"id": ..., "run": ..., "states": N or [ids]}, ...],
    "available_configurations": <unused configurations for production>, optional "available_states": <reproducible states per
    configuration>, optional "grid": {"deltat": .., "nstep": ..}}. Run paths are relative to the manifest.
    Returns (runs, ids, expected state counts, available, grid, available_states)."""
    m = read_json(path)
    rows = m.get("configurations")
    if not isinstance(rows, list) or not rows:
        raise ConfigError("manifest: 'configurations' must be a non-empty list")
    runs, ids, counts = [], [], []
    for r in rows:
        if not {"id", "run", "states"} <= set(r):
            raise ConfigError(f"manifest: every configuration needs id, run and states, got {r!r}")
        runs.append(str((Path(path).parent / r["run"]).resolve()))
        ids.append(str(r["id"]))
        counts.append(len(r["states"]) if isinstance(r["states"], list) else int(r["states"]))
    if len(set(ids)) != len(ids):
        raise ConfigError("manifest: configuration IDs repeat")
    return runs, ids, counts, m.get("available_configurations"), m.get("grid"), m.get("available_states")


def expected_run_shape(s):
    """(states, time points) every run must have, from base_input (the states of configuration.s_range, A94; nstep); None
    where not fixed."""
    conf = s.base_input.get("configuration", {})
    sr, nstep = (conf.get("s_range") if isinstance(conf, dict) else None), s.base_input.get("nstep")
    states = sr[1] - sr[0] + 1 if isinstance(sr, list) and len(sr) == 2 and all(isinstance(x, int) for x in sr) else None
    return (states if states and states > 0 else None, nstep if isinstance(nstep, int) and nstep > 0 else None)
