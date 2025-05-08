{
    # Config - Method 
    "method"        : "CCE", 
    "quantity"      : "coherence",
    
    # Config - File
    "qubitfile"     :  "/home/huijin/cal/bath/5.DiaExspins/P1/P1bath_10ppm/bath_DiaP1_10ppm_Defect",
    "gyrofile"      :  "/home/huijin/cal/bath/5.DiaExspins/P1/P1_gyro",
    "bathfile"      : ["/home/huijin/cal/bath/5.DiaExspins/P1/P1bath_10ppm/bath_DiaP1_10ppm_"],
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

,  "Defect"         : [
		{
			"dfname"	: "P1",
			"naddspin" 	: 1,
			"types" 	: ["15N"],
			"spins" 	: [0.5],   
			"gyros" 	: [-2.71261804],
			"navaax" 	: 4,
		
			"hypf" 		: [ [1, "15N",  [-114.060798,     0.000000,    0.000000,    0.000000,   -114.060798,     0.000000 ,    0.000000,    0.000000,   -159.951079]],  
			  				[2, "15N",  [-144.654779,   -17.663442,  -12.489096,  -17.663442,   -124.258790,    -7.210583 ,  -12.489096,   -7.210583,   -119.159104]], 
			  				[3, "15N",  [-144.654772,    17.663446,   12.489095,   17.663446,   -124.258798,    -7.210585 ,   12.489095,   -7.210585,   -119.159104]],
			  				[4, "15N",  [-114.060798,     0.000001,    0.000000,    0.000001,   -154.852772,     14.421165,    0.000000,   14.421165,   -119.159104]]]
		}
                       ]
}
