#!/usr/bin/env python3
"""The state and axis files of a CCEX input: configuration.bath_state_file, defect_axis_file and defect_state_file.

    python3 make_state_files.py INPUT [--force]

INPUT is a CCEX input (.json or .jsonc; //, # and ! start a comment outside strings). Its configuration block names the
files: bathfile (one or more bath files), k_range (every k; {k} in the names), s_range ([first, last]: the sampled states;
{s} in a state file name is the state number, without {s} the number is appended, as CCEX reads them), seed, path (the
folder of the block's names) and gyrofile. Its Defect.species give n_axis_states, n_onsite_spins and spins. Relative names
are read from the input file's folder. For every k, every bath file whose name in the key is not "" and (state files) every
s of s_range, the tool writes:

  bath_state_file    a line with the bath file's first line, then for every bath-file line m_s drawn uniformly from the
                     2S+1 levels of its species (S from the gyrofile, else the built-in isotope table): S - floor(u (2S+1)),
                     u from numpy default_rng([seed, k, s]).random, one draw per line
  defect_axis_file   the same first line, then for every line the axis state, uniform over 1..n_axis_states of its Defect
                     species (u from default_rng([seed, k, 0])), 0 for a spin that is not a Defect species with axis states
  defect_state_file  the same first line, then for every line the m of each on-site spin of its Defect species (the
                     species' spins, uniform over their 2I+1 levels, u from default_rng([seed, k, s, 1]): one draw per
                     on-site spin, in the order of the lines), space-separated; 0 for a line of another species

seed: configuration.seed when it is a number, its entry for k when it is a list (one per configuration of k_range), 1
when absent. Without k_range k is 0. The bath-state and axis rules are those of the automatic convergence
(example/tools/AutoConvergence), which makes its files with these functions. Files that exist are not overwritten
without --force. The rule is the tool's own: CCEX draws the states it is not given with its own random numbers.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
ISOTOPE_TABLE = HERE.parents[3] / "src" / "input" / "isotope_table.cpp" if len(HERE.parents) > 3 else Path("/nonexistent")
STATE_RULE = "m_s uniform over the 2S+1 levels, numpy default_rng([seed, k, s]), one draw per bath-file line"
AXIS_RULE = "axis uniform over 1..n_axis_states of the Defect species, numpy default_rng([seed, k, 0]), 0 elsewhere"
ONSITE_RULE = "m of each on-site spin uniform over its 2I+1 levels, numpy default_rng([seed, k, s, 1]), one draw per on-site spin"


class Refused(ValueError):
    """An input the tool cannot make files for, with the reason."""


def strip_comments(text):
    """The input without //, # and ! comments (outside strings), as CCEX reads it."""
    out, in_str, esc, i = [], False, False, 0
    while i < len(text):
        c = text[i]
        if in_str:
            out.append(c)
            if esc:
                esc = False
            elif c == "\\":
                esc = True
            elif c == '"':
                in_str = False
        elif c in "#!" or text.startswith("//", i):
            while i < len(text) and text[i] != "\n":
                i += 1
            continue
        else:
            if c == '"':
                in_str = True
            out.append(c)
        i += 1
    return "".join(out)


def read_input(path):
    return json.loads(strip_comments(Path(path).read_text()))


def get(d, key, default=None):
    """d[key] looked up without regard to case, as CCEX reads keys."""
    if isinstance(d, dict):
        for k, v in d.items():
            if k.lower() == key.lower():
                return v
    return default


def name_for(name, k=None, s=None):
    """A file name of the configuration block with {k} (with k_range) and {s} (a state number) put in; {{k}} and {{s}} are
    the literal {k} and {s}. A state file name without {s} gets the number appended."""
    out, i, found = [], 0, False
    while i < len(name):
        for key, val in (("k", k), ("s", s)):
            if val is not None and name.startswith("{{" + key + "}}", i):
                out.append("{" + key + "}"); i += len(key) + 4; break
            if val is not None and name.startswith("{" + key + "}", i):
                out.append(str(val)); i += len(key) + 2; found = found or key == "s"; break
        else:
            out.append(name[i]); i += 1
    text = "".join(out)
    return text if s is None or found else text + str(s)


def bath_lines(bathfile):
    """(the bath file's first line, the species name of every spin line)."""
    lines = Path(bathfile).read_text().splitlines()
    names = [r[3] if len(r) > 3 else "" for r in (line.split() for line in lines[1:]) if r]
    return lines[0].split()[0] if lines and lines[0].split() else str(len(names) + 1), names


def species_spins(names, gyrofile=None, table=ISOTOPE_TABLE):
    """name -> S: the gyrofile's lines (name spin gyro) first, then the built-in isotope table."""
    out = {}
    if gyrofile is not None:
        for line in Path(gyrofile).read_text().splitlines():
            f = line.split()
            if len(f) >= 3 and f[0] in names:
                out.setdefault(f[0], float(f[1]))
    missing = [n for n in names if n not in out]
    if missing and Path(table).is_file():
        for m in re.finditer(r'\{"([^"]+)",\s*([-0-9.eE+]+)f?,\s*([-0-9.eE+]+)\}', Path(table).read_text()):
            if m.group(1) in missing:
                out[m.group(1)] = float(m.group(2))
    return out


def levels(S, u):
    """m drawn from the 2S+1 levels S, S-1, ..., -S with u in [0, 1)."""
    return S - np.floor(u * (2 * S + 1))


def bath_states(names, spins, seed, k, s):
    S = np.array([spins[n] for n in names], dtype=float)
    return levels(S, np.random.default_rng([int(seed), int(k), int(s)]).random(len(names)))


def axes(names, n_axis, seed, k):
    """n_axis: species name (lower case: CCEX compares Defect names ignoring case) -> n_axis_states."""
    n = np.array([n_axis.get(x.lower(), 0) for x in names])
    u = np.random.default_rng([int(seed), int(k), 0]).random(len(names))
    return np.where(n > 0, 1 + np.floor(u * n), 0).astype(int)


def onsite_states(names, onsite, seed, k, s):
    """onsite: species name (lower case) -> the spins of its on-site nuclear spins. A list of rows, one per line."""
    rows = [onsite.get(x.lower(), []) for x in names]
    u = np.random.default_rng([int(seed), int(k), int(s), 1]).random(sum(len(r) for r in rows))
    out, j = [], 0
    for r in rows:
        out.append([float(levels(I, u[j + i])) for i, I in enumerate(r)])
        j += len(r)
    return out


def write_bath_states(path, bathfile, spins, seed, k, s):
    head, names = bath_lines(bathfile)
    Path(path).write_text(f"{head}\n" + "".join(f"{v:.5f}\n" for v in bath_states(names, spins, seed, k, s)))


def write_axes(path, bathfile, n_axis, seed, k):
    head, names = bath_lines(bathfile)
    Path(path).write_text(f"{head}\n" + "".join(f"{v}\n" for v in axes(names, n_axis, seed, k)))


def write_onsite_states(path, bathfile, onsite, seed, k, s):
    head, names = bath_lines(bathfile)
    Path(path).write_text(f"{head}\n" + "".join((" ".join(f"{v:.5f}" for v in r) if r else "0") + "\n"
                                                for r in onsite_states(names, onsite, seed, k, s)))


def _int_pair(v, key, lowest):
    if not (isinstance(v, list) and len(v) == 2 and all(isinstance(x, int) and not isinstance(x, bool) and x >= lowest for x in v)
            and v[0] <= v[1]):
        raise Refused(f"configuration.{key} must be [first, last], two integers >= {lowest} with first <= last")
    return v


def plan(inp, root):
    """Every file to write: (kind, path, bathfile, k, s, seed)."""
    conf = get(inp, "configuration")
    if not isinstance(conf, dict):
        raise Refused("the input has no configuration block")
    if get(conf, "nstate") is not None:
        raise Refused("configuration.nstate is now configuration.s_range: [1, N] (example/tools/InputKeys/convert_input_keys.py)")
    folder = get(conf, "path")
    def place(name):
        p = Path(name if folder is None or Path(name).is_absolute() else str(Path(folder) / name))
        return p if p.is_absolute() else Path(root) / p
    kr = get(conf, "k_range")
    ks = list(range(_int_pair(kr, "k_range", 0)[0], kr[1] + 1)) if kr is not None else [0]
    sr = get(conf, "s_range")
    ss = list(range(_int_pair(sr, "s_range", 1)[0], sr[1] + 1)) if sr is not None else []
    seed = get(conf, "seed", 1)
    if isinstance(seed, list):
        if len(seed) != len(ks) or not all(isinstance(x, int) and not isinstance(x, bool) for x in seed):
            raise Refused(f"configuration.seed as a list holds one integer per configuration ({len(ks)})")
        seeds = dict(zip(ks, seed))
    elif isinstance(seed, int) and not isinstance(seed, bool):
        seeds = {k: seed for k in ks}
    else:
        raise Refused("configuration.seed must be an integer or a list of one per configuration")
    baths = get(conf, "bathfile")
    if not (isinstance(baths, list) and baths and all(isinstance(b, str) for b in baths)):
        raise Refused("configuration.bathfile must be a list of file names")
    kname = (lambda n, k: name_for(n, k)) if kr is not None else (lambda n, k: n)
    out = []
    for key in ("bath_state_file", "defect_axis_file", "defect_state_file"):
        names = get(conf, key)
        if names is None:
            continue
        if not (isinstance(names, list) and len(names) == len(baths) and all(isinstance(n, str) for n in names)):
            raise Refused(f"configuration.{key} must be a list of one name per bath file ({len(baths)})")
        if key != "defect_axis_file" and not ss:
            print(f"note: configuration.{key} is not written: without configuration.s_range the run is the ensemble "
                  "(no states)", file=sys.stderr)
            continue
        for k in ks:
            for b, n in zip(baths, names):
                if not n:
                    continue
                bf = place(kname(b, k))
                for s in (ss if key != "defect_axis_file" else [None]):
                    out.append((key, place(name_for(kname(n, k), s=s)), bf, k, s, seeds[k]))
    return out


def defect_tables(inp):
    """(n_axis: name -> n_axis_states, onsite: name -> on-site spins), names in lower case."""
    species = get(get(inp, "Defect"), "species", [])
    n_axis, onsite = {}, {}
    for d in species if isinstance(species, list) else []:
        name = str(get(d, "name", "")).lower()
        n_axis[name] = int(get(d, "n_axis_states", 0) or 0)
        n = int(get(d, "n_onsite_spins", 0) or 0)
        sp = get(d, "spins", [])
        if n > 0:
            if not (isinstance(sp, list) and len(sp) == n):
                raise Refused(f"Defect species {name}: spins must list the {n} on-site spins")
            onsite[name] = [float(x) for x in sp]
    return {n: v for n, v in n_axis.items() if v > 0}, onsite


def make(input_path, force=False, table=ISOTOPE_TABLE):
    """Write every file of the input; returns the paths written."""
    input_path = Path(input_path).resolve()
    inp, root = read_input(input_path), input_path.parent
    todo = plan(inp, root)
    n_axis, onsite = defect_tables(inp)
    exists = [str(p) for _, p, *_ in todo if p.exists()]
    if exists and not force:
        raise Refused(f"{len(exists)} files exist (first {exists[0]}): --force to overwrite")
    conf = get(inp, "configuration")
    gf = get(conf, "gyrofile")
    if gf is not None:
        folder = get(conf, "path")
        gf = Path(gf if folder is None or Path(gf).is_absolute() else str(Path(folder) / gf))
        gf = gf if gf.is_absolute() else root / gf
    spins = {}
    written = []
    for key, path, bf, k, s, seed in todo:
        if not bf.is_file():
            raise Refused(f"bath file {bf} (configuration k = {k}) does not open")
        path.parent.mkdir(parents=True, exist_ok=True)
        if key == "bath_state_file":
            names = set(bath_lines(bf)[1])
            missing = names - set(spins)
            if missing:
                spins.update(species_spins(missing, gf, table))
                unknown = sorted(missing - set(spins))
                if unknown:
                    raise Refused(f"species {unknown}: neither in the gyrofile nor in the built-in isotope table ({table})")
            write_bath_states(path, bf, spins, seed, k, s)
        elif key == "defect_axis_file":
            if not n_axis:
                raise Refused("configuration.defect_axis_file: no Defect species has n_axis_states > 0")
            write_axes(path, bf, n_axis, seed, k)
        else:
            if not onsite:
                raise Refused("configuration.defect_state_file: no Defect species has n_onsite_spins > 0")
            write_onsite_states(path, bf, onsite, seed, k, s)
        written.append(path)
    return written


def main(argv):
    args = [a for a in argv if a != "--force"]
    if len(args) != 1 or args[0] in ("-h", "--help"):
        print(__doc__)
        return 0 if args and args[0] in ("-h", "--help") else 2
    try:
        written = make(args[0], force="--force" in argv)
    except (Refused, OSError, ValueError) as e:
        print(f"make_state_files: {e}", file=sys.stderr)
        return 2
    print(f"{len(written)} files written" + (f" ({written[0]} ... {written[-1]})" if written else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
