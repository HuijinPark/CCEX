{
    # Config (General)
    "method"        : "gCCE", 
    "quantity"      : "coherence",
    "qubitfile"     : "../../bath/bath_Diamond_NV_natab_defect",
    "gyrofile"      : "../../bath/gyro_13C",
    "bathfile"      : ["../../bath/bath_Diamond_NV_natab_1"],
    "order"         : 2,
    "bfield"        : 500,
    "rbath"         : 50,
    "rdip"          : 10,
    "deltat"        : 0.001,
    "nstep"         : 2000,
    "nstate"        : 1,
    "seed"          : 9407,

    # Qubit 
    "alphams"     : 1,
    "betams"      : 0,

    # Pulse
    "npulse"      : 4,
    "pulsename"   : "WAHUHA",

    # Output
    "outfile"     : "./results/CCE2_Diamond_NV_natab_1"

}
