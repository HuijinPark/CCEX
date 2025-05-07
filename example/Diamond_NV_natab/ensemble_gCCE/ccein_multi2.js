{
    # Config (General)
    "method"        : "gCCE", 
    "quantity"      : "coherence",
    "qubitfile"     : "../bath/bath_Diamond_NV_natab_defect",
    "gyrofile"      : "../bath/gyro_13C",
    "bathfile"      : ["../bath/bath_Diamond_NV_natab_1","../bath/bath_Diamond_NV_natab_1"],
    "order"         : 2,
    "bfield"        : 500,
    "rbath"         : 50,
    "rdip"          : 10,
    "deltat"        : 0.001,
    "nstep"         : 2000,
    "bathadjust"    : [[0,0,0],[0,0,50]]

    # Qubit 
    "Qubit"   : { 
                   "nqubit" : 2,
                   "qubit"  : [{   
                                    "name"   : "qubit1",
                                    "spin"   : 1.0,
                                    "gyro"   : -17608.597050,
                                    "xyz"    :  [  73.73895,	  75.72808,	  73.97033],
                                    "alphams": -1,
                                    "betams" : 0,
                                    "detuning": [0.0]
                              },
                              {   
                                    "name"   : "qubit2",
                                    "spin"   : 1.0,
                                    "gyro"   : -17608.597050,
                                    "xyz"    :  [  73.73895,	  75.72808,	  123.97033],
                                    "alphams": -1,
                                    "betams" : 0,
                                    "detuning": [0.0]
                              }
                   ],  
        			"intmap": [{
        			                "between": ["qubit1", "qubit1"],
        			                "tensor" : [[-960000, 0, 0],[0, -960000, 0 ],[0, 0, 1920000]]
        			           },
                               {
        			                "between": ["qubit1", "qubit1"],
        			                "tensor" : [[-960000, 0, 0],[0, -960000, 0 ],[0, 0, 1920000]]
        			           }
                    ]
                },   

    "alphams"     : 1,
    "betams"      : 0,

    # Pulse
    "npulse"      : 1,

    # Output
    "outfile"     : "./results/CCE2_Diamond_NV_natab_1"

}
