{
    # Config (General)
    "method"        : "gCCE", 
    "quantity"      : "dm",
    "gyrofile"      : "./gyro_13C",
    "bathfile"      : ["./bath"],
    "order"         : 1,
    "bfield"        : 500,
    "rbath"         : 500,
    "rdip"          : 100,
    "deltat"        : 0.001,
    "nstep"         : 2000,
    "nstate"        : 1,

    "Qubit": {
        "nqubit": 1,
        "qubit": [
        {
            "name": "qubit1",
            "spin": 1.0,
            "gyro": -17608.597050,
            "xyz": [0.0,0.0,0.0],
            "alpha": [0.0],
            "beta": [1.0],
            "detuning": [0.0]
        }]
    },

    # Qubit 
    #"alphams"     : 1,
        #"betams"      : 0,

    # Pulse
    "npulse"      : 1,

    # Output
    "outfile"     : "./results/CCE2_Diamond_NV_natab_1",
    "savemode"      : "allfull"
}
