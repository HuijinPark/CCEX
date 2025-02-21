{
    # Config - Method 
    "method"        : "CCE", 
    "quantity"      : "coherence",
    
    # Config - File
    "qubitfile"     :  "/home/huijin/tutorial/Diamond_NV_P1Bath/bath/14N_1ppm/bath_DiaP1_1ppm_Defect",
    "gyrofile"      :  "/home/huijin/tutorial/Diamond_NV_P1Bath/bath/DiaP1_gyro",
    "bathfile"      : ["/home/huijin/tutorial/Diamond_NV_P1Bath/bath/14N_1ppm/bath_DiaP1_1ppm_"],
  #"avaaxfile"     :  "Random",
  #"statefile"     :  "Random",
  #"exstatefile"   :  "Random",

    # Config - General
    "order"         : 2,
    "bfield"        : 500,
    "rbath"         : 850,
    "rdip"          : 620,
    "deltat"        : 0.01,
    "nstep"         : 300,  

    # Qubit 
    "alphams"     : -1,
    "betams"      : 0,

    # Cluster
    "addsubclus"  : 0, # true / false
    "nk"          : [0,0,0], # default : [0] (all) or [0,0,30,40,50]

    # Pulse
    "npulse"      : 1,

    # Config - Spin tensor
    "hf_readmode"   : 0, 
#  "hf_tensorfile" : "0",
#  "hf_cutoff"     : 0
#  "qd_readmode"   : 0,

    # Output
    "savemode"    : "normal"

}
