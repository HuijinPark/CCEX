{
    # Config (General)
    "method"        : "CCE",
    "quantity"      : "coherence",
    "qubitfile"     : "../bath/bath_Diamond_NV_natab_defect",
    "gyrofile"      : "../bath/gyro_13C",
    "bathfile"      : ["../bath/bath_Diamond_NV_natab_cubic_1"],
    "order"         : 2,
    "bfield"        : 500,
    "rbath"         : 50,
    "rdip"          : 10,
    "deltat"        : 0.001,
    "nstep"         : 2000,

    # Qubit : ms = +1 / 0 of the NV, one pi pulse (Hahn echo)
    "qalphams"       : 1,
    "qbetams"        : 0,
    "npulse"        : 1,

    # Coordinate frame : the bath file is written in the cubic crystal frame, where the NV
    # axis is [111]. Rotate positions (about the reference qubit) so that qubit_axis lands
    # on the z axis, along which "bfield" is applied. "enabled": false or no block at all
    # leaves the bath as written.
    "coordinate_frame_rotation" : {
        "enabled"    : true,
        "bath_axis"  : [0, 0, 1],
        "qubit_axis" : [1, 1, 1]
    },

    "outfile"       : "./results/CCE2_cubic_rotated_1"
}
