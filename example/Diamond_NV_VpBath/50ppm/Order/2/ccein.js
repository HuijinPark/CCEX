{
    # Config - Method 
    "method"        : "CCE", 
    "quantity"      : "coherence",
    
    # Config - File
    "qubitfile"     :  "/home/huijin/git/CCEX/example/Diamond_NV_VpBath/bath/50ppm/bath_DiaVp_50ppm_Defect",
    "gyrofile"      :  "/home/huijin/git/CCEX/example/Diamond_NV_VpBath/bath/DiaVp_gyro",
    "bathfile"      : ["/home/huijin/git/CCEX/example/Diamond_NV_VpBath/bath/50ppm/bath_DiaVp_50ppm_"],
  #"avaaxfile"     :  "Random",
  #"statefile"     :  "Random",
  #"exstatefile"   :  "Random",

    # Config - General
    "order"         : 2,
    "bfield"        : 500,
    "rbath"         : 240,
    "rdip"          : 220,
    "deltat"        : 0.0002,
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
