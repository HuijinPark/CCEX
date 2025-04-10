{
    "method" : "gcce",
    "quantity" : "coherence",

    "gyrofile" : "./gyro",
    "bathfile" : ["./bath_1"],

    "hf_readmode" : 3,
    "hf_tensorfile" : "./Atensor",

    # Qubit #(unit : kHz),  D = 2.88GHz
    "Qubit"   : {
                   "nqubit" : 1,
                   "qubit"  : [
                   {
                       "name"   : "qubit1",
                       "spin"   : 1.0,
                       "gyro"   : -17608.597050,
                       "xyz"    :  [  523.336, 524.936, 521.003],
                       "alphams": -1,
                       "betams" : 0,
                       "detuning": [0.0]
                   }],


                    "intmap": [
                    {
                        "between": ["qubit1", "qubit1"],
                        "tensor": [
                            [-960000, 0, 0],
                            [0, -960000, 0 ],
                            [0, 0, 1920000 ]
                        ]
                    }]
    },

    "order"    : 1,
    "npulse"   : 1,
    #"pulsename": "HahnEcho",
    #"sequence" : [
    #              ["1/4", "Y", 180], 
    #              ["2/4", "Y", 180],
    #              ["1/8", "Y",  180],
    #              ["1/8", "I",   0]
    #             ],
    "bfield" : 448.84,
    "deltat": 0.00001,
    "nstep": 50, 

    "rbath": 10,
    "rdip": 10,
    
    "outfile" : "./CPMG_result"

}
