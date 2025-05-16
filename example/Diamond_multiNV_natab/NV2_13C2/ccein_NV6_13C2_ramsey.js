{
    "method" : "gCCE",
    "quantity" : "dm",
    "order" : 2,
    "nstep" : 300,
    "deltat" : 0.0000001,

    "Qubit" : {
                "nqubit" : 2,
                "qubit" : [
                    {"name": "q1", "xyz":[ 75.663000, 74.262800, 74.134800], "spin": 1.0, "gyro": -17608.597050, "alphams": 1, "betams": 0},
                    {"name": "q2", "xyz":[ 75.663000, 74.262800,124.134800], "spin": 1.0, "gyro": -17608.597050, "alphams": 1, "betams": 0}
                ],
                "intmap" : [
                    {"between" : ["q1","q1"], "tensor":[[-960000,0,0],[0,-960000,0],[0,0,1920000]]},
                    {"between" : ["q2","q2"], "tensor":[[-960000,0,0],[0,-960000,0],[0,0,1920000]]}
                ]
    },

    "bathfile"  : "./bath",
    "statefile" : "./state",
    "gyrofile"  : "./gyro_13C",
    "rbath"     : 10000,
    "rdip"      : 50000,
    "nstate"    : 1,
    "bfield"    : 500,

    "npulse"      : 0,
    
    "outfile"     : "./results/CCE2_Diamond_6NV_13C2_1"


}
