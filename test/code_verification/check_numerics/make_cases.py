"""make_cases.py: the cases and their inputs, cut from example/tutorial. Run from the repository root; it
rewrites inputs/ and cases/ (not reference/). Sizes are chosen so that every case runs in seconds; the
state and axis files are written here with a fixed pattern so that two binaries read the same files."""
import os, re, json, shutil, glob
T = 'example/tutorial'; R = 'test/code_verification/check_numerics'
def lenient(t):
    t = re.sub(r'[#!][^\n]*', '', t); t = re.sub(r',\s*([}\]])', r'\1', t); return json.loads(t)
FILE_KEYS = ('bathfile', 'gyrofile', 'qubitfile', 'hf_tensorfile', 'qd_tensorfile', 'qd_tensorfile_woqubit', 'avaaxfile', 'statefile', 'exstatefile')
# name: (source input, system folder under inputs/, overrides, deterministic files wanted)
C = {}
def case(name, src, sys_, ov, det=()): C[name] = dict(src=f'{T}/{src}', sys=sys_, ov=ov, det=det)
nat = dict(rbath=30, nstep=400, seed=19)
case('natab_cce2_hahn',        '02_Diamond_NV_natab/ensemble_CCE/ccein.json',           'natab', nat)
case('natab_cce2_ramsey',      '02_Diamond_NV_natab/ensemble_CCE/ccein.json',           'natab', dict(nat, npulse=0))
case('natab_cce2_cpmg8',       '02_Diamond_NV_natab/ensemble_CCE/ccein.json',           'natab', dict(nat, npulse=8))
case('natab_cce2_single',      '02_Diamond_NV_natab/single-sample_CCE/ccein.json',      'natab', dict(nat, nstate=2), ('state',))
# The two gCCE2 cases were cut from ensemble_gCCE/ccein_single.json and ccein_multi2.json, whose qubitfile won over an
# unread Qubit block (one NV in both). Those inputs changed on 2026-10-13 (a qubitfile with a Qubit block is refused;
# multi2 is a real two-qubit input, single has the ZFS), so the cases are cut from ensemble_CCE's input with the same
# values: the same inputs as recorded, without the Qubit block. natab_gcce2_bathadjust (natab_gcce2_two_qubits until
# 2026-10-13) is the one NV with the bath file given twice, the second copy shifted 50 A along z.
gc = dict(method='gCCE', rbath=20, nstep=300, seed=19)
case('natab_gcce2_hahn',       '02_Diamond_NV_natab/ensemble_CCE/ccein.json',           'natab', gc)
case('natab_gcce2_bathadjust', '02_Diamond_NV_natab/ensemble_CCE/ccein.json',           'natab',
     dict(gc, bathfile=['../bath/bath_Diamond_NV_natab_1', '../bath/bath_Diamond_NV_natab_1'], bathadjust=[[0, 0, 0], [0, 0, 50]]))
case('natab_dsj2',             '02_Diamond_NV_natab/ensemble_Dsj/ccein.json',           'natab', nat)
case('natab_knight_shift',     '02_Diamond_NV_natab/knight_shift/ccein_wi_knight.json', 'natab', nat)
case('natab_manual_hahn',      '06_Pulse_Manual/HahnEcho_natab_ensemble_origin/ccein.json', 'natab_pulse', nat)
p1 = dict(rbath=500, nstep=200, seed=19)
case('p1_cce2_ensemble',       '04_Diamond_NV_P1Bath/SpinEcho_CCE2_1ppm/ccein.json',    'p1_1ppm', dict(p1, nstate=0), ('avaax',))
case('p1_cce2_single',         '04_Diamond_NV_P1Bath/SpinEcho_CCE2_1ppm/ccein.json',    'p1_1ppm', dict(p1, nstate=2), ('avaax', 'state', 'exstate'))
case('p1_cce2e2n',             '04_Diamond_NV_P1Bath/SpinEcho_CCE2e2n_1ppm/ccein.json', 'p1_1ppm', dict(rbath=400, nstep=200, seed=19, nstate=1), ('avaax', 'state', 'exstate'))
case('p1_cce4',                '04_Diamond_NV_P1Bath/SpinEcho_CCE4_1ppm/ccein.json',    'p1_1ppm', dict(rbath=500, nstep=200, seed=19, nstate=1), ('avaax', 'state', 'exstate'))
case('p1_gcce2_zfs',           '04_Diamond_NV_P1Bath/SpinEcho_gCCE2_10ppm_ZFS/ccein_with_ZFS.json', 'p1_10ppm', dict(rbath=200, seed=19, nstate=1), ('avaax', 'state', 'exstate'))
case('p1_pcce2_sk2',           '04_Diamond_NV_P1Bath/SpinEcho_pCCE2_1ppm/ccein_sK2.json', 'p1_1ppm', dict(p1, nstate=1), ('avaax', 'state', 'exstate'))
case('p1_gcce2_population',    '04_Diamond_NV_P1Bath/Population_gCCE2_10ppm/ccein.json', 'p1_10ppm', dict(rbath=250, nstep=500, seed=19, nstate=1), ('avaax', 'state', 'exstate'))
case('p1_cpmg8_manual',        '06_Pulse_Manual/CPMG_P1_1ppm/ccein_origin.json',        'p1_1ppm', dict(p1, nstate=1), ('avaax', 'state', 'exstate'))
case('p1_manual_hahn',         '06_Pulse_Manual/SpinEcho_P1_1ppm/ccein_origin.json',    'p1_1ppm', dict(p1, nstate=1), ('avaax', 'state', 'exstate'))
case('vp_cce2',                '05_Diamond_NV_VpBath/SpinEcho_CCE2/ccein_10ppm.json',   'vp_10ppm', dict(rbath=200, nstep=200, seed=19, nstate=1), ('avaax', 'state', 'exstate'))
case('hbn_cce2_hf3_qd2',       '07_hBN_VB_natab/ensemble_CCE/ccein_hf3_qd2.json',       'hbn', dict(seed=19))
case('hbn_cce2_hf3_qd0',       '07_hBN_VB_natab/ensemble_CCE/ccein_hf3_qd0.json',       'hbn', dict(seed=19))
case('hbn_hfmedi',             '07_hBN_VB_natab/hyperfine_mediated_term/ccein_wi_hfmedi.json', 'hbn', dict(seed=19))
case('hbn_gcce2_dm',           '07_hBN_VB_natab/single-sample_gCCE/ccein.json',         'hbn', dict(seed=19, nstate=1), ('state',))
case('hbn_gcce2_oh',           '07_hBN_VB_natab/single-sample_gCCE_DPS/ccein_wi_OH.json', 'hbn', dict(seed=19, nstate=1), ('state',))
case('arbitrary_gcce1_dm',     '01_Arbitrary_system/ccein.json',                        'arbitrary', dict(seed=19, nstate=1, nstep=200), ('state',))

shutil.rmtree(f'{R}/inputs', ignore_errors=True); shutil.rmtree(f'{R}/cases', ignore_errors=True); os.makedirs(f'{R}/inputs'); os.makedirs(f'{R}/cases')
def bath_lines(path):
    """(count, [isotope names]) of a bath file: first line is the count, the rest 'x y z name'."""
    rows = [l.split() for l in open(path).read().splitlines()[1:] if l.strip()]
    return len(rows), [r[3] if len(r) > 3 else '' for r in rows]
def spins_from_gyro(path):
    s = {}
    for l in open(path):
        p = l.split()
        if len(p) >= 3:
            try: s[p[0]] = float(p[1])
            except ValueError: pass
    return s
def det_files(sysdir, bath, gyro, defect, want):
    """avaax_ (one file), state_1..2 and exstate_1..2 with a fixed pattern; every file is indexed by the bath file's lines."""
    n, names = bath_lines(bath); spins = spins_from_gyro(gyro) if gyro else {}
    made = {}
    if 'avaax' in want:
        navaax = int(defect.get('navaax', 0)) if defect else 0
        if navaax > 0:
            p = f'{sysdir}/avaax_nav{navaax}'; open(p, 'w').write(f'{n}\n' + ''.join(f'{(i * 7) % navaax + 1}\n' for i in range(n))); made['avaaxfile'] = f'avaax_nav{navaax}'
    if 'state' in want:
        for k in (1, 2):
            vals = []
            for i, nm in enumerate(names):
                S = spins.get(nm, 0.5)
                # a valid sublevel of this spin: the first two files alternate the sign
                if S == int(S): vals.append(0.0)
                else: vals.append(0.5 if (i + k) % 2 else -0.5)
            open(f'{sysdir}/state_{k}', 'w').write(f'{n}\n' + ''.join(f'{v:.5f}\n' for v in vals))
        made['statefile'] = 'state_'
    if 'exstate' in want and defect and int(defect.get('naddspin', 0)) > 0:
        for k in (1, 2):
            open(f'{sysdir}/exstate_{k}', 'w').write(f'{n}\n' + ''.join(f'{float(((i + k) % 3) - 1):.5f}\n' for i in range(n)))
        made['exstatefile'] = 'exstate_'
    return made
summary = []
for name, c in C.items():
    src_dir = os.path.dirname(c['src']); cfg = lenient(open(c['src'], errors='replace').read())
    sysdir = f'{R}/inputs/{c["sys"]}'; os.makedirs(sysdir, exist_ok=True)
    cfg.pop('outfile', None)
    if isinstance(cfg.get('bathfile'), str): cfg['bathfile'] = [cfg['bathfile']]     # a bare string is read as no bath at all
    for k, v in c['ov'].items(): cfg[k] = v
    # every file the input names is copied once into inputs/<system>/ and the input points there
    def bring(ref, key=''):
        srcp = os.path.normpath(os.path.join(src_dir, ref))
        if os.path.isfile(srcp):
            dst = f'{sysdir}/{os.path.basename(srcp)}'
            if not os.path.exists(dst): shutil.copy(srcp, dst)
            return f'../../inputs/{c["sys"]}/{os.path.basename(srcp)}'
        # a prefix: a bath prefix (Runner style, "bath_DiaP1_1ppm_") means configuration 1; a state prefix
        # ("state") means the files state1, state2 the calculation will index
        if os.path.isfile(srcp + '1') and key == 'bathfile':
            dst = f'{sysdir}/{os.path.basename(srcp)}1'
            if not os.path.exists(dst): shutil.copy(srcp + '1', dst)
            return f'../../inputs/{c["sys"]}/{os.path.basename(srcp)}1'
        got = [g for g in (srcp + '1', srcp + '2') if os.path.isfile(g)]
        for g in got:
            dst = f'{sysdir}/{os.path.basename(g)}'
            if not os.path.exists(dst): shutil.copy(g, dst)
        return f'../../inputs/{c["sys"]}/{os.path.basename(srcp)}' if got else None
    for key in FILE_KEYS:
        if key not in cfg: continue
        v = cfg[key]
        if isinstance(v, list): cfg[key] = [bring(x, key) or x for x in v]
        elif isinstance(v, str) and v not in ('Random',):
            nv = bring(v, key)
            if nv: cfg[key] = nv
            else: del cfg[key]                       # 'Random' or a file that never existed: the deterministic file below replaces it
    defect = (cfg.get('Defect') or [None])[0] if isinstance(cfg.get('Defect'), list) else cfg.get('Defect')
    bath = os.path.join(R, 'cases', name, cfg['bathfile'][0] if isinstance(cfg['bathfile'], list) else cfg['bathfile'])
    bath = os.path.normpath(os.path.join(f'{R}/cases/{name}', cfg['bathfile'][0] if isinstance(cfg['bathfile'], list) else cfg['bathfile']))
    gyro = os.path.normpath(os.path.join(f'{R}/cases/{name}', cfg['gyrofile'])) if 'gyrofile' in cfg else None
    if c['det']:
        made = det_files(sysdir, bath, gyro, defect, c['det'])
        for k, f in made.items(): cfg[k] = f'../../inputs/{c["sys"]}/{f}'
    os.makedirs(f'{R}/cases/{name}', exist_ok=True)
    json.dump(cfg, open(f'{R}/cases/{name}/ccein.json', 'w'), indent=2)
    summary.append((name, cfg.get('method'), cfg.get('quantity'), cfg.get('nstate', 0), cfg.get('npulse'), cfg.get('rbath'), cfg.get('nstep'), sorted(k for k in FILE_KEYS if k in cfg)))
for s in summary: print(f'{s[0]:24} {str(s[1]):5} {str(s[2]):9} nstate={s[3]} npulse={s[4]} rbath={s[5]} nstep={s[6]} files={",".join(s[7])}')
print('inputs:', sum(len(f) for _, _, f in os.walk(f'{R}/inputs')), 'files,', end=' '); os.system(f'du -sh {R}/inputs | cut -f1')
