#!/usr/bin/env python3
"""
Axis (2) : write the new-form inputs next to the legacy input they must match.

Each variant is a copy of the system's legacy ccein with ONE thing changed, so a
failure names the key that caused it. The manifest that comes out says, per
variant, which baseline it is compared against and how:

    byte  the two inputs must mean exactly the same number to the parser, so the
          raw output has to be byte-identical
    num   a unit conversion or a constructed tensor sits between them; the two
          reach the same value by different arithmetic, so the comparison reports
          the largest deviation instead of demanding byte identity

Nothing here stores an expected coherence. Every check is legacy-vs-new inside
one binary, and axis (1) is what ties the legacy side back to v1.1.0.
"""
import json, sys, os, copy

sysname, work = sys.argv[1], sys.argv[2]
legacy = json.load(open(f"{work}/legacy.json"))
rows = []

def emit(name, d, baseline="legacy", mode="byte", desc=""):
    d = copy.deepcopy(d)
    d["outfile"] = f"./{name}"
    json.dump(d, open(f"{work}/{name}.json", "w"), indent=4)
    rows.append((name, baseline, mode, desc))

# ---------------------------------------------------------------- rotation keys
ROT_OFF = {"enabled": False, "bath_axis": [0.0,0.0,1.0], "qubit_axis": [1.0,1.0,1.0]}
ROT_ID  = {"enabled": True,  "bath_axis": [0.0,0.0,1.0], "qubit_axis": [0.0,0.0,1.0]}
ROT_111 = {"enabled": True,  "bath_axis": [0.0,0.0,1.0], "qubit_axis": [1.0,1.0,1.0]}

def with_rot(section, key="coordinate_frame_rotation", **extra):
    d = copy.deepcopy(legacy); d[key] = dict(section); d[key].update(extra); return d

# The baseline every variant is measured against is the legacy run itself, which
# axis (1) already showed is byte-identical to v1.1.0.
emit("base", legacy, desc="legacy input, unchanged")

emit("rot_disabled", with_rot(ROT_OFF), "base", "byte",
     'coordinate_frame_rotation.enabled=false must be a no-op')

# The bare identity variant is only valid where no further key is REQUIRED. S4 reads
# a DFT tensor file, so hf/qd_tensor_frame has to be stated; S5 has two qubits, so
# reference_qubit has to be. Both refusals are the point of those keys, and each of
# those systems emits its own identity case below with them supplied.
if sysname == "S4":
    pass
elif sysname == "S5":
    pass
else:
    emit("rot_identity", with_rot(ROT_ID), "base", "byte",
         'identity rotation (qubit_axis=[0,0,1]) must be a no-op')

# ---------------------------------------------------------------- per system
if sysname == "S1a":
    # The deprecated name has to mean exactly what the new one means.
    emit("rot111_newname", with_rot(ROT_111), "base", "num",
         '[111] rotation : moves the bath, so it is NOT expected to match the baseline')
    emit("rot111_alias", with_rot(ROT_111, key="bath_coordinate_rotation"),
         "rot111_newname", "byte",
         'deprecated bath_coordinate_rotation == coordinate_frame_rotation')

elif sysname == "S3":
    df = legacy["Defect"][0]
    hypf_rows = df["hypf"]; efg_rows = df["efg"]
    gyro = df["gyros"][0];  eq = df["eqs"][0]
    navaax = df["navaax"]

    def mk(**defect_changes):
        d = copy.deepcopy(legacy)
        nd = d["Defect"][0]
        # object forms all require the frame to be explicit; "qubit" is the legacy
        # default, so stating it must not move a number
        nd["coordinate_frame"] = "qubit"
        nd.update(defect_changes)
        return d

    emit("df_frame_only", mk(), "base", "byte",
         'Defect.coordinate_frame="qubit" is the legacy default')

    # canonical units : the scale factor is exactly 1, so byte identity is required
    emit("df_gyro_canon", mk(gyros={"values":[gyro], "unit":"radkHz/G"}),
         "df_frame_only", "byte", 'gyros {values,unit=radkHz/G} == legacy array')
    emit("df_eq_canon",   mk(eqs={"values":[eq], "unit":"1e-30 m^2"}),
         "df_frame_only", "byte", 'eqs {values,unit=1e-30 m^2} == legacy array')
    emit("df_hypf_canon", mk(hypf={"values":hypf_rows, "unit":"MHz"}),
         "df_frame_only", "byte", 'hypf {values,unit=MHz} == legacy array')
    emit("df_efg_canon",  mk(efg={"values":efg_rows, "unit":"Hartree/Bohr^2"}),
         "df_frame_only", "byte", 'efg {values,unit=Hartree/Bohr^2} == legacy array')

    # a different unit for the same physical quantity : conversion, so compare numerically
    S_GYRO = 2.0*3.141592653589793*1.0e-1          # MHz/T   -> radkHz/G
    S_EQ   = 1.0e-1                                 # mbarn   -> 1e-30 m^2
    S_HF   = 1.0e3                                  # GHz     -> MHz
    hartree_to_ev = 27.211386; bohr_to_m = 0.5291772e-10
    S_EFG  = 1.0e20*bohr_to_m*bohr_to_m/hartree_to_ev   # V/A^2 -> Hartree/Bohr^2

    def scaled(rows_, s):
        return [[r[0], r[1], [v/s for v in r[2]]] for r in rows_]

    emit("df_gyro_mhzt", mk(gyros={"values":[gyro/S_GYRO], "unit":"MHz/T"}),
         "df_frame_only", "num", 'gyros in MHz/T == the same gyro in radkHz/G')
    emit("df_eq_mbarn",  mk(eqs={"values":[eq/S_EQ], "unit":"mbarn"}),
         "df_frame_only", "num", 'eqs in mbarn == the same eQ in 1e-30 m^2')
    emit("df_hypf_ghz",  mk(hypf={"values":scaled(hypf_rows,S_HF), "unit":"GHz"}),
         "df_frame_only", "num", 'hypf in GHz == the same tensors in MHz')
    emit("df_efg_va2",   mk(efg={"values":scaled(efg_rows,S_EFG), "unit":"V/angstrom^2"}),
         "df_frame_only", "num", 'efg in V/angstrom^2 == the same tensors in a.u.')

    # detuning : absent from this fixture, so add it in both forms at once
    det_rows = [[i, "e", 0.35*i] for i in range(1, navaax+1)]   # arbitrary, fixed
    d_leg = copy.deepcopy(legacy); d_leg["Defect"][0]["detuning"] = det_rows
    emit("df_det_legacy", d_leg, "base", "num",
         'detuning added (legacy array) : changes the physics, its own baseline')
    emit("df_det_object", mk(detuning={"values":det_rows, "unit":"MHz"}),
         "df_det_legacy", "byte", 'detuning {values,unit=MHz} == legacy indexed array')
    emit("df_det_ghz", mk(detuning={"values":[[r[0],r[1],r[2]/1.0e3] for r in det_rows],
                                    "unit":"GHz"}),
         "df_det_legacy", "num", 'detuning in GHz == the same values in MHz')

    # defect_axis_reference is documented as changing no number; require exactly that
    d = with_rot(ROT_ID); d["defect_axis_reference"] = "qubit_axis"
    emit("defaxis_ref", d, "rot_identity", "byte",
         'defect_axis_reference="qubit_axis" is validated and logged, nothing else')

elif sysname == "S4":
    # With the rotation ON but R = I, "bath" and "qubit" must both reduce to the
    # unrotated answer. This pins the PARSING and the no-op path; the [111] case for
    # these same keys is what test/rotation/run.sh [13,14,15] and its A/B section do.
    emit("tf_bath_identity",
         with_rot(ROT_ID, hf_tensor_frame="bath", qd_tensor_frame="bath"),
         "base", "byte", 'identity rotation with hf/qd_tensor_frame="bath"')
    emit("tf_qubit_identity",
         with_rot(ROT_ID, hf_tensor_frame="qubit", qd_tensor_frame="qubit"),
         "base", "byte", 'identity rotation with hf/qd_tensor_frame="qubit"')

elif sysname == "S5":
    def with_intmap_frame(frame):
        d = copy.deepcopy(legacy)
        for e in d["Qubit"]["intmap"]:
            e["tensor_frame"] = frame
        return d
    emit("intmap_qubit_frame", with_intmap_frame("qubit"), "base", "byte",
         'intmap tensor_frame="qubit", rotation off : declaring the frame changes nothing')

    # The identity case needs both keys nqubit>1 makes mandatory: reference_qubit, and
    # a tensor_frame on every explicit intmap entry. Refusing it without them is the
    # behaviour the [10] refusal tests cover; here it is supplied and must be a no-op.
    d = with_intmap_frame("bath")
    d["coordinate_frame_rotation"] = dict(ROT_ID, reference_qubit="q1",
                                          qubit_position_frame="bath")
    emit("rot_identity", d, "base", "byte",
         'identity rotation (reference_qubit=q1, intmap frames given) must be a no-op')
    emit("intmap_bath_frame", with_intmap_frame("bath"), "base", "byte",
         'intmap tensor_frame="bath", rotation off : same')

    # reference_qubit moves the rotation ORIGIN only. R is the same in both, so the two
    # coordinate sets differ by the rigid translation (I-R)(r0-r1) and the coherence
    # must not move. The first run would break if the key were ignored.
    for q in ("q1", "q2"):
        d = with_intmap_frame("bath")
        d["coordinate_frame_rotation"] = dict(ROT_111,
                                              reference_qubit=q,
                                              qubit_position_frame="bath")
        emit(f"refq_{q}", d, "refq_q1", "num",
             f'[111] rotation about {q}')

with open(f"{work}/manifest.tsv", "w") as fh:
    for r in rows:
        fh.write("\t".join(r) + "\n")
print(f"    {len(rows)} variants")
