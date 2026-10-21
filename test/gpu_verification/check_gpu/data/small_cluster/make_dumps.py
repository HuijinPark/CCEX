#!/usr/bin/env python3
"""The CPU-dumped clusters that test_small_cluster_cuda compares the batched small-cluster GPU kernel against (A26).

  python3 make_dumps.py BIN [OUT]                 the ten fixture files (single/ensemble x Hahn/CPMG/Ramsey x exact/recursive)
  python3 make_dumps.py --step0 CLUSTERS.BIN [OUT] p1_production_cr60_single_hahn_exact.bin

BIN is a CPU build (bin/main.out, any preset; 2026-10-12: wsl-gnu, commit 1d6550cc0). The fixture is
test/code_verification/check_input_forms/inputs/Diamond_NV_P1_10ppm_cce without its Defect block, so that a P1 is the
electron alone: order 1 gives dimension 2, order 2 dimension 4. Each run, one rank, writes the first 150 clusters of
dimension 2 and 4 (55 + 95) with CCEX_SMALL_CLUSTER_DUMP (small_cluster_batch.h). CLUSTERS.BIN is the A26 step-0 dump of
the production P1 CR60 sigma4 cCCE2 input (rank 0, state 1; idm4-gpu job 674,
~/tmp/ccex-a26-20261012/prod/CR060_sigma04_cCCE2_dump/clusters.bin): its first 32 clusters of dimension 2 and 96 of
dimension 4 are rewritten in the format of the dump hook (single state, Hahn, exact, 151 points).
p1_pcce2_spins_single_hahn_recursive.bin (A30, 2026-10-13) was written by the hook directly: check_rdip_ladder's
p1_pcce2_spins base input (check_numerics p1_pcce2_sk2 with pcce.connections spins), 1 rank, CCEX_SMALL_CLUSTER_DUMP_N=150.
"""
import json, os, shutil, struct, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = HERE
CASES = {
    'single_hahn_exact': dict(nstate=1, pulse={'count': 1}),
    'single_hahn_recursive': dict(nstate=1, pulse={'count': 1}, evolution_update='recursive'),
    'ensemble_hahn_exact': dict(nstate=0, pulse={'count': 1}),
    'ensemble_hahn_recursive': dict(nstate=0, pulse={'count': 1}, evolution_update='recursive'),
    'single_cpmg4_exact': dict(nstate=1, pulse={'count': 4}),
    'single_cpmg3_recursive': dict(nstate=1, pulse={'count': 3}, evolution_update='recursive'),
    'ensemble_cpmg2_exact': dict(nstate=0, pulse={'count': 2}),
    'ensemble_cpmg3_recursive': dict(nstate=0, pulse={'count': 3}, evolution_update='recursive'),
    'single_ramsey_exact': dict(nstate=1, pulse={'count': 0}),
    'ensemble_ramsey_recursive': dict(nstate=0, pulse={'count': 0}, evolution_update='recursive'),
}


def fixture(binary, out):
    src = os.path.join(HERE, '../../../../code_verification/check_input_forms/inputs/Diamond_NV_P1_10ppm_cce/inputfiles')
    work = tempfile.mkdtemp(prefix='small_cluster_')
    shutil.copytree(src, os.path.join(work, 'inputfiles'))
    base = json.load(open(os.path.join(src, 'ccein_nstate1.json')))
    del base['Defect']
    for name, change in CASES.items():
        c = dict(base); c.update(change); c['output'] = {'file': './' + name}
        c['configuration'] = dict(base['configuration'], s_range=[1, c.pop('nstate')])   # configuration.s_range (A64, A94)
        if 'evolution_update' in c:                                                 # the execution block's (A65)
            c['execution'] = dict(base.get('execution', {}), evolution_update=c.pop('evolution_update'))
        json.dump(c, open(os.path.join(work, name + '.json'), 'w'), indent=1)
        target = os.path.join(out, name + '.bin')
        if os.path.exists(target): os.remove(target)
        env = dict(os.environ, CCEX_SMALL_CLUSTER_DUMP=target, CCEX_SMALL_CLUSTER_DUMP_N='150')
        subprocess.run(['mpirun', '-n', '1', binary, '-f', name + '.json'], cwd=work, env=env, check=True,
                       stdout=subprocess.DEVNULL)
        print(name, os.path.getsize(target), 'bytes')
    shutil.rmtree(work)


def step0(clusters, out, n2=32, n4=96):
    f = open(clusters, 'rb').read(); o = 0; taken = {2: 0, 4: 0}; want = {2: n2, 4: n4}
    with open(os.path.join(out, 'p1_production_cr60_single_hahn_exact.bin'), 'wb') as w:
        while o < len(f):
            d, ns = struct.unpack_from('ii', f, o); dt, = struct.unpack_from('d', f, o + 8); o += 16
            body = f[o:o + 16 * (2 * d * d + d)]; o += 16 * (2 * d * d + d)
            cpu = f[o:o + 16 * ns]; o += 16 * ns
            if taken[d] >= want[d]: continue
            taken[d] += 1
            # d, ensemble, recursive, intervals, distinct, nstep; deltat; fraction[distinct]; same[intervals]
            w.write(struct.pack('6i', d, 0, 0, 2, 1, ns) + struct.pack('d', dt) + struct.pack('d', 0.5) + struct.pack('2i', 0, 0))
            w.write(body + cpu)
    print('p1_production_cr60_single_hahn_exact', taken)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--step0':
        step0(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else OUT)
    elif len(sys.argv) > 1:
        fixture(os.path.abspath(sys.argv[1]), sys.argv[2] if len(sys.argv) > 2 else OUT)
    else:
        sys.exit(__doc__)
