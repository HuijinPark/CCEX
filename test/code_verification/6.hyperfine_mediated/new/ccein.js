{
    # Config (General)
    "method"        : "CCE", 
    "quantity"      : "coherence",
    "qubitfile"     : "../../CCE_Reprod/Bath_Data/4.hexagonal_data/defect",
    "gyrofile"      : "../../CCE_Reprod/Bath_Data/2.ensemble_data/Gyro_h_BN",
        #"bathfile"      : ["../bath/flat/bath_hBN_VB_flat_1"],
    "order"         : 2,
    "bfield"        : 30000,
    "rbath"         : 20,
    "rdip"          : 10,
    "deltat"        : 0.002,
    "nstep"         : 100,
    "nstate"        : 0,


    "qzfs"          : [[-1080000,0,0],[-1180000,0,0],[0,0,2320000]],
    "hfmedi"        : true,

    # Qubit 
    "alphams"     : 1,
    "betams"      : 0,

    # Pulse
    "npulse"      : 1,

    # Bath HF
    "hf_readmode" : 3,
    "hf_tensorfile" : "../../CCE_Reprod/Bath_Data/4.hexagonal_data/Afile_vertex",
    "hf_cutoff" : 0,
    "hf_ignore_oor" : 0,

    # Bath QD
    "qd_readmode" : 2,
    "qd_tensorfile" : "../../CCE_Reprod/Bath_Data/4.hexagonal_data/Qfile_vertex",
    "qd_tensorfile_woqubit" : "../../CCE_Reprod/Bath_Data/4.hexagonal_data/Qfile_vertex",

    # Output
    "outfile"     : "./gCCE2_hBN_VB_natab_1"

}
