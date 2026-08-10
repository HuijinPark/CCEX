{
   "method"        : "gCCE",
   "quantity"      : "coherence",
   "seed"          : 3528,
   "gyrofile"      : "GYROFILE",

   "Qubit" : {
      "nqubit" : 1,
      "qubit"  : [
         { "name" : "qubit1", "spin" : 1.0, "gyro" : -17608.597050,
            "xyz" : [0.0, 0.0, 0.0], "alphams" : -1, "betams" : 0, "detuning": [0.0] } ],
      "intmap" : [
         { "between" : ["qubit1", "qubit1"],
            "tensor"  : [ [-960000.0, 0, 0], [0, -960000.0, 0], [0, 0, 1920000.0] ] } ],
      "overhaus" : true
   },

   "order"         : 1,
   "bfield"        : 180,
   "rbath"         : 200,
   "rdip"          : 200,
   "deltat"        : 2.5e-04,
   "nstep"         : 41,
   "alphams"       : -1,
   "betams"        : 0,
   "addsubclus"    : 0,
   "nk"            : [0,0],
   "npulse"        : 1,
   "hf_readmode"   : 0,
   "savemode"      : "normal",
   "bnpulse"             : 1,
   "bpulse_defect"       : "P1",
   "bpulse_energy_shift" : [79.975539],
   "bspin"               : 0.5,
   "balphams"            : 0.5,
   "bbetams"             : -0.5,
   "bsequence"           : [["0.50", "X", 180.0], ["0.50", "I", 0.0]],

   "Defect" : [
      { "dfname"   : "P1",
         "apprx"    : true,
         "naddspin" : 0,
         "navaax"   : 8,
         "detuning" : [
            [  1, "e",    79.975539],
            [  2, "e",    59.579552],
            [  3, "e",    59.579552],
            [  4, "e",    59.579552],
            [  5, "e",   -79.975539],
            [  6, "e",   -59.579552],
            [  7, "e",   -59.579552],
            [  8, "e",   -59.579552] ]
      }
   ]
}
