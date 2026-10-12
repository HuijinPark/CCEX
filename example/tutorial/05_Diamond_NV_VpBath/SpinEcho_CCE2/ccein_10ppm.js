{
    # Config - Method 
    "method"        : "CCE", 
    "quantity"      : "coherence",
    
    # Config - File
    "qubitfile"     :  "../bath/10ppm/bath_DiaVp_10ppm_Defect",
    "gyrofile"      :  "../bath/DiaVp_gyro",
    "bathfile"      : ["../bath/10ppm/bath_DiaVp_10ppm_"],
  #"avaaxfile"     :  "Random",
  #"statefile"     :  "Random",
  #"exstatefile"   :  "Random",

    # Config - General
    "order"         : 2,
    "bfield"        : 500,
    "rbath"         : 370,
    "rdip"          : 330,
    "deltat"        : 0.00003,
    "nstep"         : 300,  

    # Qubit 
    "qalphams"     : -1,
    "qbetams"      : 0,

    # Cluster
    "addsubclus"  : 0, # true / false
    "nk"          : [0,0,0], # default : [0] (all) or [0,0,30,40,50]

    # Pulse
    "npulse"      : 1,

    # Config - Spin tensor
    "hf_readmode"   : 0, 
#  "hf_tensorfile" : "0",
#  "qd_readmode"   : 0,

    # Output
    "savemode"    : "normal"

}
