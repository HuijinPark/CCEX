{
    # Config (General)
    "method"        : "CCE", 
    "quantity"      : "coherence",
    "qubitfile"     : "../bath/flat/bath_hBN_VB_flat_Defect",
    "gyrofile"      : "../bath/hBN_natab_gyro",
    "bathfile"      : ["../bath/flat/bath_hBN_VB_flat_1"],
    "order"         : 2,
    "bfield"        : 30000,
    "rbath"         : 20,
    "rdip"          : 10,
    "deltat"        : 0.003,
    "nstep"         : 50,
    "nstate"        : 0,

    # Qubit 
    "alphams"     : 1,
    "betams"      : 0,

    # Pulse
    "npulse"      : 1,

    # Bath HF
    "hf_readmode"   : 3,
    "hf_tensorfile" : "./../bath/flat/hBN_VB_flat_ATensor_HSE_CellPBE_tag",
    "hf_cutoff"     : 0,
    "hf_ignore_oor" : 0,
    "DefectTotSpin" : 1.0,
    "CorrTotSpin"   : 0.0,

    # Bath QD
    "qd_readmode" : 0, #0 (X), exp, dft(2 widefect) dft(wi+wo defect), dft(wi+wo defect)
    "qd_tensorfile" : "./../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE",
    "qd_tensorfile_woqubit" : "./../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE",
    "qd_cellpara" : [0.0, 0.0, 0.0],

    # Output
    "outfile"     : "./gCCE2_hBN_VB_natab_1"

    # Bath
    #"Bath" : {
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] }
    #}

}
