"""The figures of docs/user/Examples.md, from the results the tutorials write.

Run the tutorial calculations first (each folder's run.sh; see the table in
docs/user/Examples.md for which ones and with what NP/NCONF/NSTATE), then

    python3 docs/images/tutorials/make_figures.py        # needs numpy, matplotlib

from the repository root. A figure whose results are missing is skipped and
named. Coherence figures show |L(t)|; density-matrix figures show the diagonal
of rho(t) and the largest off-diagonal element. The time column of a CCEX output is in ms; the unit of each figure
is the one the tutorial's own post-processing uses."""
import glob, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'example', 'analysis', 'Analyzer', 'CoherenceAnalyzer.v1'))
from coherence import read_coherence
T = os.path.join(ROOT, 'example', 'tutorial')
OUT = os.path.dirname(os.path.abspath(__file__))
SCALE = {'ms': 1.0, 'us': 1e3, 'ns': 1e6}
UNIT = {'ms': 'ms', 'us': 'µs', 'ns': 'ns'}


def read_dm(path):
    """A density-matrix output: blocks of rows, the first row of a block starting
    with the time. Returns (t, rho) with rho of shape (nstep, d, d)."""
    times, blocks, current = [], [], []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            if not parts[0].endswith('j'):
                if current:
                    blocks.append(current)
                times.append(float(parts[0])); current = []
                parts = parts[1:]
            current.extend(complex(p) for p in parts)
    if current:
        blocks.append(current)
    d = int(round(len(blocks[0]) ** 0.5))
    return np.array(times), np.array(blocks).reshape(len(blocks), d, d)


FIGURES = [
    ('01_arbitrary_dm', 'dm', 'ms', '01 gCCE1 density matrix of a hand-written system',
     [('rho_ii(t)', '01_Arbitrary_system/results/CCE2_Diamond_NV_natab_1_0_wiDiv')]),
    ('02_ensemble_CCE', 'L', 'ms', '02 NV in natural-abundance diamond, CCE2 Hahn echo',
     [('CCE2, 1 configuration', '02_Diamond_NV_natab/ensemble_CCE/results/CCE2_Diamond_NV_natab_EnsAvg_conv1')]),
    ('02_knight_shift', 'L', 'us', '02 Knight shift in a tilted field',
     [('without knight', '02_Diamond_NV_natab/knight_shift/results/CCE2_Diamond_NV_natab_wo_knight_1_wiDiv'),
      ('with knight', '02_Diamond_NV_natab/knight_shift/results/CCE2_Diamond_NV_natab_wi_knight_1_wiDiv')]),
    ('03_NV2_dm', 'dm', 'ns', '03 two NV centres with a 13C pair, gCCE density matrix, Hahn echo',
     [('rho_ii(t)', '03_Diamond_multiNV_natab/NV2_13C2/results/CCE2_Diamond_NV_natab_NV6_13C2_hahnecho_1_state1_wiDiv')]),
    ('04_P1_orders', 'L', 'ms', '04 NV in a 1 ppm P1 bath, spin echo by order',
     [('CCE2 (50 states)', '04_Diamond_NV_P1Bath/SpinEcho_CCE2_1ppm/results/CCE2_DiaP1_1ppm_Order_2_EnsAvg_conv50'),
      ('CCE4 (50 states)', '04_Diamond_NV_P1Bath/SpinEcho_CCE4_1ppm/results/CCE4_DiaP1_1ppm_Order_4_EnsAvg_conv50'),
      ('CCE2e2n (5 states)', '04_Diamond_NV_P1Bath/SpinEcho_CCE2e2n_1ppm/results/CCE2e2n_DiaP1_1ppm_Order_2_EnsAvg_conv5')]),
    ('04_gCCE2_ZFS', 'L', 'ms', '04 gCCE2 at zero field, 10 ppm, with and without the zero-field splitting',
     [('with ZFS', '04_Diamond_NV_P1Bath/SpinEcho_gCCE2_10ppm_ZFS/results/gCCE2_DiaP1_10ppm_B0_0G_with_ZFS_EnsAvg_conv1'),
      ('without ZFS', '04_Diamond_NV_P1Bath/SpinEcho_gCCE2_10ppm_ZFS/results/gCCE2_DiaP1_10ppm_B0_0G_without_ZFS_EnsAvg_conv1')]),
    ('04_population', 'dm', 'ms', '04 gCCE2 population of the NV in a 10 ppm P1 bath',
     [('rho_ii(t)', '04_Diamond_NV_P1Bath/Population_gCCE2_10ppm/results/gCCE2_DiaP1_10ppm_1_wiDiv')]),
    ('05_Vp_ppm', 'L', 'ms', '05 NV in a Vp bath, CCE2 spin echo by concentration',
     [('10 ppm', '05_Diamond_NV_VpBath/SpinEcho_CCE2/results/CCE2_DiaVp_10ppm_Order_2_EnsAvg_conv1'),
      ('12 ppm', '05_Diamond_NV_VpBath/SpinEcho_CCE2/results/CCE2_DiaVp_12ppm_Order_2_EnsAvg_conv1'),
      ('50 ppm', '05_Diamond_NV_VpBath/SpinEcho_CCE2/results/CCE2_DiaVp_50ppm_Order_2_EnsAvg_conv1')]),
    ('06_CPMG', 'L', 'ms', '06 CPMG with eight pulses on the 1 ppm P1 bath, three spellings, 10 sampled states',
     [('origin', '06_Pulse_Manual/CPMG_P1_1ppm/results/CCE2_DiaP1_1ppm_Order_2_origin_EnsAvg_conv10'),
      ('test', '06_Pulse_Manual/CPMG_P1_1ppm/results/CCE2_DiaP1_1ppm_Order_2_test_EnsAvg_conv10'),
      ('wahuha-style', '06_Pulse_Manual/CPMG_P1_1ppm/results/CCE2_DiaP1_1ppm_Order_2_wahuha_EnsAvg_conv10')]),
    ('06_hahnecho_manual', 'L', 'ms', '06 Hahn echo written as a sequence, and WAHUHA (gCCE2)',
     [('Hahn echo, sequence, ensemble', '06_Pulse_Manual/HahnEcho_natab_ensemble_origin/results/CCE2_Diamond_NV_natab_EnsAvg_conv1'),
      ('WAHUHA, gCCE2, one state', '06_Pulse_Manual/WAHUHA_natab_single_gcce/results/CCE2_Diamond_NV_natab_EnsAvg_conv1')]),
    ('07_hBN_tensors', 'L', 'us', '07 V_B in h-BN, CCE2 with the DFT hyperfine (hf) and quadrupole (qd) tensors',
     [('hf 0, qd 2', '07_hBN_VB_natab/ensemble_CCE/results/CCE2_hBN_VB_hf0_qd2_1_wiDiv'),
      ('hf 3, qd 0', '07_hBN_VB_natab/ensemble_CCE/results/CCE2_hBN_VB_hf3_qd0_1_wiDiv'),
      ('hf 3, qd 2', '07_hBN_VB_natab/ensemble_CCE/results/CCE2_hBN_VB_hf3_qd2_1_wiDiv')]),
    ('07_hfmedi', 'L', 'us', '07 the hyperfine-mediated term',
     [('without hfmedi', '07_hBN_VB_natab/hyperfine_mediated_term/results/CCE2_hBN_VB_wo_hfmedi_1_wiDiv'),
      ('with hfmedi', '07_hBN_VB_natab/hyperfine_mediated_term/results/CCE2_hBN_VB_wi_hfmedi_1_wiDiv')]),
    ('07_DPS_OH', 'L', 'ns', '07 gCCE2 with and without the OH term (DPS)',
     [('without OH', '07_hBN_VB_natab/single-sample_gCCE_DPS/results/CCE2_hBN_VB_wo_OH_1_state1_wiDiv'),
      ('with OH', '07_hBN_VB_natab/single-sample_gCCE_DPS/results/CCE2_hBN_VB_wi_OH_1_state1_wiDiv')]),
    ('08_bundle', 'L', 'us', '08 CCE2 with ordinary clusters and with bundles (bundle_max_dim 8)',
     [('standard clusters', '08_Bundle_clustering/CCE2_bundle/results/CCE2_standard_1_wiDiv'),
      ('bundles', '08_Bundle_clustering/CCE2_bundle/results/CCE2_bundle_1_wiDiv')]),
    ('09_backend', 'L', 'us', '09 CCE5 of one 1024-dimensional cluster, per backend',
     [('backend cpu', '09_GPU_backend/CCE5_P1/results/CCE5_P1_cpu_1_wiDiv'),
      ('backend auto', '09_GPU_backend/CCE5_P1/results/CCE5_P1_auto_1_wiDiv')]),
    ('10_frame', 'L', 'us', '10 the bath in the cubic frame, with and without the rotation',
     [('aligned bath', '10_Coordinate_rotation/HahnEcho_CCE2/results/CCE2_aligned_1_wiDiv'),
      ('cubic frame, rotation on', '10_Coordinate_rotation/HahnEcho_CCE2/results/CCE2_cubic_rotated_1_wiDiv'),
      ('cubic frame, rotation off', '10_Coordinate_rotation/HahnEcho_CCE2/results/CCE2_cubic_norotation_1_wiDiv')]),
]


def resolve(rel):
    hits = sorted(glob.glob(os.path.join(T, rel)))
    return hits[0] if hits else None


def main():
    made, skipped = [], []
    for name, kind, unit, title, curves in FIGURES:
        paths = [(label, resolve(rel)) for label, rel in curves]
        if any(p is None for _, p in paths):
            skipped.append((name, [rel for (label, rel), (_, p) in zip(curves, paths) if p is None]))
            continue
        fig, ax = plt.subplots(figsize=(6.4, 4.0))
        for label, path in paths:
            if kind == 'L':
                t, L = read_coherence(path, time_scale=SCALE[unit])
                ax.plot(t, np.abs(L), label=label, lw=1.4)
            else:
                t, rho = read_dm(path)
                for i in range(rho.shape[1]):
                    ax.plot(t * SCALE[unit], rho[:, i, i].real, label=f'rho_{i}{i}', lw=1.2)
                # the coherence sits off the diagonal: the element largest at t = 0
                off = np.abs(rho[0]) * (1 - np.eye(rho.shape[1]))
                if off.max() > 0:
                    i, j = np.unravel_index(off.argmax(), off.shape)
                    ax.plot(t * SCALE[unit], np.abs(rho[:, i, j]), '--', label=f'|rho_{i}{j}|', lw=1.2)
                ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel(f'time ({UNIT[unit]})')
        ax.set_ylabel('|L(t)|' if kind == 'L' else 'population')
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, name + '.png'), dpi=150)
        plt.close(fig)
        made.append(name)
    print('made', len(made), ':', ' '.join(made))
    for name, missing in skipped:
        print('skipped', name, '- missing', missing)


if __name__ == '__main__':
    main()
