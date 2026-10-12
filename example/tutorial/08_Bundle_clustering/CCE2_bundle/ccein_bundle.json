{
    # Config (General)
    "method"        : "CCE",
    "quantity"      : "coherence",
    "qubitfile"     : "../bath/bath_Diamond_NV_natab_defect",
    "gyrofile"      : "../bath/gyro_13C",
    "bathfile"      : ["../bath/bath_Diamond_NV_natab_1"],
    "order"         : 2,
    "bfield"        : 500,
    "rbath"         : 30,
    "rdip"          : 10,
    "deltat"        : 0.001,
    "nstep"         : 2000,

    # Clustering : "cluster_mode": "adaptive" ties strongly coupled spins into bundles (docs/clustering-terminology.md);
    #              order counts bundles per cluster, bundle_max_dim caps one bundle's Hilbert dimension
    #              (8 = up to three 13C). Bundles of one spin each are ordinary CCE2. With rdip present
    #              the bundles are connected by distance; without it weak connections are dropped after a check.
    "cluster_mode"    : "adaptive",
    "bundle_max_dim"  : 8,

    # Qubit : ms = +1 / 0 of the NV, one pi pulse (Hahn echo)
    "qalphams"       : 1,
    "qbetams"        : 0,
    "npulse"        : 1,

    "outfile"       : "./results/CCE2_bundle_1"
}
