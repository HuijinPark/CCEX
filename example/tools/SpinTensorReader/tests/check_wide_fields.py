#!/usr/bin/env python3
"""Values that fill the 10-character field (2026-10-13, (internal notes) "텐서 도구").

The writers put "{:10.5f}" fields between tabs. CCEX's tensor reader split a line at spaces only, so a value of 10
or more characters (-123.45678, 1234.56789, -12345.60000) ran into the field before it ("-0.09147\\t-141.13231\\t")
and the row was dropped (fewer than 13 numbers) or read shifted. Since then the writers put a space in front of every
value and the reader splits at any blank.

For the VASP and the QE writer, this rewrites the h-BN A and Q files of check_numerics with the writers' own line
formats (taken from the module source, so no numpy is needed), puts large values into the first tensor row, an eQ
line and a quadrupole etc line, and runs CCEX -v on the hbn_cce2_hf3_qd2 case: every number CCEX prints back for the
A tensor, the Q tensor, the eQ and the Q etc lines must be the number written. Each writer's file is checked twice:
as written now, and with the format before the fix (no space in front of the values), which the reader must read too.

Usage:  check_wide_fields.py <main.out> <mpirun> <work dir>      (exit 0 = pass)
"""
import ast, json, os, re, subprocess, sys

BIN, MPIRUN, WORK = sys.argv[1:4]
TOOL = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.abspath(os.path.join(TOOL, "../../.."))
CASE = os.path.join(REPO, "test/code_verification/check_numerics/cases/hbn_cce2_hf3_qd2")
INP = os.path.join(REPO, "test/code_verification/check_numerics/inputs/hbn")
A_SRC = os.path.join(INP, "hBN_VB_flat_ATensor_HSE_CellPBE_v3")
Q_SRC = os.path.join(INP, "hBN_VB_flat_QTensor_HSE_CellPBE")


def formats(module):
    """The line-format functions of a writer module, without importing the package."""
    src = open(os.path.join(TOOL, "spHamilReader", module)).read()
    names = {"divline", "Afilehead", "Afileinfo", "Atensorform", "Qfilehead", "Qfileinfo", "Qtensorform"}
    tree = ast.parse(src)
    tree.body = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    ns = {}
    exec(compile(tree, module, "exec"), ns)
    return ns


def legacy(s):
    """The format before the fix: no space in front of the values."""
    return s.replace(" {:10.5f}", "{:10.5f}").replace("\t {:>10f}", "\t{:>10f}")


def parse(path, ncol):
    """Tag/provenance lines, g-factor or eQ pairs, etc rows, min/max, and the tensor rows of a tensor file."""
    pre, const, etc, rng, rows = [], [], [], [], []
    for line in open(path).read().split("\n"):
        t = line.split()
        if not t or set(line.strip()) == {"-"} or t[0] == "AtomName":
            continue
        if t[1:2] in (["g-factor___:"], ["eQ___:"]):
            const.append([t[0], float(t[2])])
        elif t[1:2] == ["etc____:"]:
            etc.append([t[0]] + [float(x) for x in t[5:]])
        elif t[0] in ("MinDif[A]:", "MaxDif[A]:"):
            rng.append([float(x) for x in t[1:]])
        elif len(t) == ncol + 1:
            rows.append([t[0]] + [float(x) for x in t[1:]])
        else:
            pre.append(line + "\n")
    return pre, const, etc, rng, rows


def write_a(f, path, old):
    pre, g, etc, rng, rows = parse(A_SRC, 13)
    rows[0][4], rows[0][5], rows[0][13] = 1234.56789, -12345.6, -123.45678        # iso, Axx, Azz
    d = f["Afileinfo"]()
    fix = legacy if old else (lambda s: s)
    with open(path, "w") as o:
        o.writelines(pre)
        o.write(f["divline"]() + f["Afilehead"]() + f["divline"]())
        o.writelines(fix(d["gline"]).format(*x) for x in g)
        o.writelines(fix(d["etcline"]).format(x[0]) for x in etc)
        o.write(fix(d["minline"]).format(*rng[0]) + fix(d["maxline"]).format(*rng[1]) + f["divline"]())
        o.writelines(fix(f["Atensorform"]()).format(*r) for r in rows)
    return g, etc, rows


def write_q(f, path, old):
    pre, eq, etc, rng, rows = parse(Q_SRC, 12)
    rows[0][4], rows[0][12] = -1234.56789, 12345.6                               # EFGxx, EFGzz
    eq[0][1] = 12345.678
    etc[0][1], etc[0][9] = -1234.5, 100.25
    d = f["Qfileinfo"]()
    fix = legacy if old else (lambda s: s)
    with open(path, "w") as o:
        o.writelines(pre)
        o.write(f["divline"]() + f["Qfilehead"]() + f["divline"]())
        o.writelines(fix(d["eQline"]).format(*x) for x in eq)
        o.writelines(fix(d["etcline"]).format(*x) for x in etc)
        o.write(fix(d["minline"]).format(*rng[0]) + fix(d["maxline"]).format(*rng[1]) + f["divline"]())
        o.writelines(fix(f["Qtensorform"]()).format(*r) for r in rows)
    return eq, etc, rows


def block(log, title):
    """The lines after a title CCEX -v prints, up to the blank line."""
    lines = log.split("\n")
    i = next(k for k, l in enumerate(lines) if title in l) + 1
    out = []
    while i < len(lines) and lines[i].strip():
        out.append(lines[i]); i += 1
    return out


def nums(line):
    return [float(x) for x in re.findall(r"[-+]?\d+\.\d+", line)]


def check(tag, log, written, a_or_q):
    """Compare what CCEX printed with what was written (to the 5 decimals of the file)."""
    bad = []
    const, etc, rows = written
    if a_or_q == "A":
        got_rows = [nums(l)[0:13] for l in block(log, "[A-tensor] (MHz)")[1:]]
        got_const = [nums(l) for l in block(log, "[g-factor] (no-unit)")]
        got_etc = [nums(l.split(":", 1)[1]) for l in block(log, "[A etc] (MHz)")[1:]]
        want_etc = [[0.0] for _ in etc]
    else:
        got_rows = [nums(l)[0:12] for l in block(log, "[Q-tensor wi defect]")[1:]]
        got_const = [nums(l) for l in block(log, "[eQ] (Q/milibarn*10^-1)")]
        got_etc = [nums(l.split(":", 1)[1]) for l in block(log, "[Q etc] (MHz)")[1:]]
        want_etc = [x[1:] for x in etc]
    want_rows = [[round(v, 5) for v in r[1:]] for r in rows]
    if len(got_rows) != len(want_rows):
        bad.append(f"{len(got_rows)} tensor rows read, {len(want_rows)} written")
    for k, (g, w) in enumerate(zip(got_rows, want_rows)):
        if len(g) != len(w) or max(abs(a - b) for a, b in zip(g, w)) > 1e-6:
            bad.append(f"row {k + 1}: read {g}, written {w}")
    for (name, w), g in zip(const, got_const):
        if not g or abs(g[-1] - w) > 1e-6:
            bad.append(f"{name} constant: read {g}, written {w}")
    for k, (g, w) in enumerate(zip(got_etc, want_etc)):
        if len(g) != len(w) or max(abs(a - b) for a, b in zip(g, w)) > 1e-6:
            bad.append(f"etc line {k + 1}: read {g}, written {w}")
    return bad


fail = 0
case = json.load(open(os.path.join(CASE, "ccein.json")))
for module in ("vasp/mkAQtensorfile_vasp.py", "qe/mkAQtensorfile_qe.py"):
    f = formats(module)
    for old in (False, True):
        tag = f"{module.split('/')[0]} {'format before the fix' if old else 'format now'}"
        d = os.path.join(WORK, module.split("/")[0] + ("_old" if old else "_new"))
        os.makedirs(d, exist_ok=True)
        wa = write_a(f, os.path.join(d, "Afile"), old)
        wq = write_q(f, os.path.join(d, "Qfile"), old)
        j = dict(case, hf_tensorfile=os.path.join(d, "Afile"), qd_tensorfile=os.path.join(d, "Qfile"),
                 qd_tensorfile_woqubit=os.path.join(d, "Qfile"), nstep=2)
        json.dump(j, open(os.path.join(d, "ccein.json"), "w"), indent=1)
        r = subprocess.run([MPIRUN, "-n", "1", BIN, "-v", "-f", os.path.join(d, "ccein.json"), "-o", os.path.join(d, "out")],
                           cwd=CASE, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        open(os.path.join(d, "log"), "w").write(r.stdout)
        bad = [f"CCEX exited {r.returncode}"] if r.returncode else []
        try:
            bad += check(tag, r.stdout, wa, "A") + check(tag, r.stdout, wq, "Q")
        except StopIteration:
            bad.append("a block of the -v output is missing")
        if bad:
            fail += 1
            print(f"  FAIL  wide fields, {tag} (see {d}/log):")
            for b in bad[:6]:
                print("        " + b)
        else:
            print(f"  PASS  wide fields, {tag}: A (iso 1234.56789, Axx -12345.6), Q (EFG -1234.56789, 12345.6), "
                  f"eQ 12345.678 and the Q etc line are read as written")
sys.exit(1 if fail else 0)
