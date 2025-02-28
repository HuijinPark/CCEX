{
    # Config - Method 
    "method"        : "CCE", 
    "quantity"      : "coherence",
    
    # Config - File
    "qubitfile"     :  "/home/huijin/cal/bath/5.DiaExspins/Vp/Vpbath_10ppm/bath_DiaVp_10ppm_Defect",
    "gyrofile"      :  "/home/huijin/cal/bath/5.DiaExspins/Vp/Vp_gyro",
    "bathfile"      : ["/home/huijin/cal/bath/5.DiaExspins/Vp/Vpbath_10ppm/bath_DiaVp_10ppm_"],
  #"avaaxfile"     :  "Random",
  #"statefile"     :  "Random",
  #"exstatefile"   :  "Random",

    # Config - General
    "order"         : 2,
    "bfield"        : 292,
    "rbath"         : 270,
    "rdip"          : 230,
    "deltat"        : 0.001,
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
