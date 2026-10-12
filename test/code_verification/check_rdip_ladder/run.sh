#!/usr/bin/env bash
#
# rdip_ladder ((internal notes)): one order-2 run (cCCE, gCCE or pCCE) that also writes the product of
# every rung r_k < rdip must give, rung by rung, the curve of a separate run at rdip = r_k, and its main output must stay the
# run without the ladder.
# ---------------------------------------------------------------------------
# Cases (from check_numerics; rungs 200, 300, 450 below rdip 620 for P1 1 ppm, 100, 200, 300 below 330 for P1 10 ppm ZFS,
# 4, 6, 8 below 10 for 13C, 3, 5, 6.5 below 8 for h-BN): cCCE p1_cce2_single (2 states) and natab_cce2_hahn (ensemble); gCCE
# on 13C, h-BN and P1 10 ppm by its three paths - ensemble (nstate 0, evolution matrix), vector (states) and matrix (states,
# evolution matrix); pCCE2 (sK 2) on 13C, h-BN (rbath 8) and P1 1 ppm with pcce_connections centers and spins
# (improvements.md A23). Criterion (registered): |difference| <= 1e-9 on wiDiv and noDiv at every point (2 ranks; the ladder
# run re-adds zero-power single units, so the rank split and the order of products can differ from the separate runs; on a
# curve that diverges to |L| ~ 1e9 and beyond that shows as a last-bit difference above 1e-9, as between 1 and 2 ranks of one
# run - (internal notes)). Then inputs the ladder must refuse. About 7 minutes.
#
# Usage:  bash test/code_verification/check_rdip_ladder/run.sh      (bin/main.out; CCEX_BIN, CCEX_MPI as elsewhere)
#         CCEX_QUICK=1 (CI): at most 30 time points and 20 pCCE partition trials per case - the same cases, rungs and refusals.
#         CCEX_JOBS=n: n runs at a time (default 1); CI uses 2 on its 4-core runner.
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../../.." && pwd)"
# One run at a time per folder: a second run would delete the first one's work files (D29; review 20261008-003400 R3).
if [ "${CCEX_RUN_LOCK:-}" != "$HERE" ]; then
    exec 9>"$HERE/.run.lock"
    if command -v flock >/dev/null 2>&1 && ! flock -n 9; then
        echo "ERROR: another check is running in $HERE (it would delete that run's work files); wait for it or use another checkout" >&2; exit 1
    fi
    CCEX_RUN_LOCK="$HERE"; export CCEX_RUN_LOCK
fi
BIN="${CCEX_BIN:-$REPO/bin/main.out}"; MPIRUN="${CCEX_MPI:-mpirun}"
CASES="$REPO/test/code_verification/check_numerics/cases"
WORK="$HERE/_work"; rm -rf "${WORK:?}"; mkdir -p "$WORK"
python3 - "$CASES" "$WORK" <<'PY'
import json, os, sys
cases, work = sys.argv[1:]
P1, P1ZFS, C13, HBN = [200, 300, 450], [100, 200, 300], [4, 6, 8], [3, 5, 6.5]
ENS = dict(nstate=0, statefile=None)                           # gCCE ensemble: the default evolution is then matrix
PCCE = dict(method='pCCE', sK=2)
# name: (check_numerics case, changes, rungs)
plan = {'p1_cce2_single': ('p1_cce2_single', {}, P1), 'natab_cce2_hahn': ('natab_cce2_hahn', {}, C13),
        'natab_gcce2_ensemble': ('natab_gcce2_hahn', {}, C13),
        'natab_gcce2_vector': ('natab_gcce2_hahn', dict(nstate=2), C13),
        'natab_gcce2_matrix': ('natab_gcce2_hahn', dict(nstate=2, evolution='matrix'), C13),
        'natab_pcce2_centers': ('natab_cce2_hahn', dict(PCCE, pcce_connections='centers'), C13),
        'natab_pcce2_spins': ('natab_cce2_hahn', dict(PCCE, pcce_connections='spins'), C13),
        'hbn_gcce2_ensemble': ('hbn_gcce2_oh', ENS, HBN),
        'hbn_gcce2_vector': ('hbn_gcce2_oh', {}, HBN),
        'hbn_gcce2_matrix': ('hbn_gcce2_oh', dict(evolution='matrix'), HBN),
        'hbn_pcce2_centers': ('hbn_hfmedi', dict(PCCE, pcce_connections='centers', rbath=8), HBN),
        'hbn_pcce2_spins': ('hbn_hfmedi', dict(PCCE, pcce_connections='spins', rbath=8), HBN),
        'p1_gcce2_ensemble': ('p1_gcce2_zfs', ENS, P1ZFS),
        'p1_gcce2_vector': ('p1_gcce2_zfs', {}, P1ZFS),
        'p1_gcce2_matrix': ('p1_gcce2_zfs', dict(evolution='matrix'), P1ZFS),
        'p1_pcce2_centers': ('p1_pcce2_sk2', dict(pcce_connections='centers'), P1),
        'p1_pcce2_spins': ('p1_pcce2_sk2', dict(pcce_connections='spins'), P1)}
runs = []
for name, (case, change, rungs) in plan.items():
    c = json.load(open(f'{cases}/{case}/ccein.json'))
    for k, v in list(c.items()):
        if isinstance(v, str) and v.startswith('../'): c[k] = os.path.normpath(f'{cases}/{case}/{v}')
        if k == 'bathfile': c[k] = [os.path.normpath(f'{cases}/{case}/{x}') for x in v]
    for k, v in change.items():
        if v is None: c.pop(k, None)
        else: c[k] = v
    if os.environ.get('CCEX_QUICK'):     # CI: what is checked (a rung = the separate run) holds at any size; the same cases, smaller
        c['nstep'] = min(int(c.get('nstep', 30)), 30)
        if str(c.get('method', '')).lower() == 'pcce': c['max_trial'] = min(int(c.get('max_trial', 20)), 20)
    d = f'{work}/{name}'; os.makedirs(d)
    jobs = {'base': dict(c), 'ladder': dict(c, rdip_ladder=rungs)}
    jobs.update({f'sep{r:g}': dict(c, rdip=r) for r in rungs})
    for job, cfg in jobs.items():
        cfg['outfile'] = f'{d}/{job}'
        json.dump(cfg, open(f'{d}/{job}.json', 'w'), indent=1); runs.append(f'{d} {job}')
    # refusals: order 3, rung >= rdip, decreasing rungs, savemode all (cCCE); quantity dm (gCCE); pCCE order 3
    bad = {}
    if case == name:
        bad = {'bad_order3': dict(c, order=3, nk=[0, 0, 0, 0], rdip_ladder=rungs[:1]), 'bad_above': dict(c, rdip_ladder=[c['rdip']]),
               'bad_order': dict(c, rdip_ladder=rungs[::-1]), 'bad_saveall': dict(c, savemode='all', rdip_ladder=rungs[:1])}
    elif name == 'natab_gcce2_vector':
        bad = {'bad_dm': dict(c, quantity='dm', rdip_ladder=rungs[:1])}
    elif name == 'natab_pcce2_centers':
        bad = {'bad_order3': dict(c, order=3, nk=[0, 0, 0, 0], rdip_ladder=rungs[:1])}
    for job, cfg in bad.items():
        cfg['outfile'] = f'{d}/{job}'; json.dump(cfg, open(f'{d}/{job}.json', 'w'), indent=1)
open(f'{work}/runs', 'w').write('\n'.join(runs) + '\n')
PY
fail=0
# CCEX_JOBS runs (each on 2 ranks) at a time; the runs are independent, so this changes the time only. With Open MPI and
# more than one at a time, set OMPI_MCA_hwloc_base_binding_policy=none, or every mpirun binds its ranks to cores 0 and 1.
run_one() {
    (cd "$1" && timeout 600 "$MPIRUN" -n 2 "$BIN" -f "$2.json" > "$2.log" 2>&1 < /dev/null) || { echo "  FAIL  $1/$2 run: see $1/$2.log"; return 1; }
}
export -f run_one; export MPIRUN BIN
xargs -P "${CCEX_JOBS:-1}" -L 1 bash -c 'run_one "$0" "$1"' < "$WORK/runs" || fail=1
for f in "$WORK"/*/bad_*.json; do
    d="$(dirname "$f")"; b="$(basename "$f" .json)"
    if (cd "$d" && timeout 300 "$MPIRUN" -n 2 "$BIN" -f "$b.json" > "$b.log" 2>&1 < /dev/null); then
        echo "  FAIL  $(basename "$d") $b was accepted"; fail=1
    elif ! grep -q "Error: rdip_ladder" "$d/$b.log"; then
        echo "  FAIL  $(basename "$d") $b failed without an rdip_ladder message: see $d/$b.log"; fail=1
    else echo "  PASS  $(basename "$d") $b refused"; fi
done
python3 - "$WORK" <<'PY' || fail=1
import glob, os, sys
work = sys.argv[1]; bad = 0
def curve(f):
    return [complex(l.split(None, 1)[1].replace(' ', '')) for l in open(f) if l.strip()]
for d in sorted(glob.glob(f'{work}/*/')):
    case = os.path.basename(d.rstrip('/'))
    pairs = []
    for f in sorted(glob.glob(f'{d}/base_*Div')):
        pairs.append((f, f.replace('/base_', '/ladder_'), 'main output'))
    for f in sorted(glob.glob(f'{d}/sep*_*Div')):
        r = os.path.basename(f).split('_')[0][3:]
        pairs.append((f, f.replace(f'/sep{r}_', f'/ladder_rdip{r}_'), f'rung {r}'))
    if not pairs:
        print(f'  FAIL  {case}: no output'); bad = 1
    worst = {}
    for ref, got, what in pairs:
        if not os.path.exists(got):
            print(f'  FAIL  {case}: {os.path.basename(got)} missing'); bad = 1; continue
        a, b = curve(ref), curve(got)
        diff = max(abs(x - y) for x, y in zip(a, b)) if len(a) == len(b) else float('inf')
        worst[what] = max(worst.get(what, 0), diff)
    for what, diff in worst.items():
        ok = diff <= 1e-9
        bad |= not ok
        print(f'  {"PASS" if ok else "FAIL"}  {case} {what}: max |difference| {diff:.2e} over wiDiv and noDiv, all states')
sys.exit(bad)
PY
[ $fail -eq 0 ] && echo "  rdip_ladder: all checks passed" || echo "  rdip_ladder: FAILED"
exit $fail
