#!/bin/sh

variable_name="PULSE"
variable_values="1 4 12 64 128"

data_directory="./data/"
data_filename_format="CCE2_DiaP1_1ppm_Pulse_PULSE_Order_2_conv400"

result="./data/CCE2_DiaP1_1ppm_Pulse_PULSE_Order_2_FittingParams"

space=5
fit_initp=1

###############################
# Single plot with ftting
python CoherenceAnalyzer.v1/run.py -d "${data_directory}" -fi "${data_filename_format}" -vn "${variable_name}" -v ${variable_values} -uc "ns" -ntn -fo "${result}" -ylim 0 1 -fit -pl_T -pl_p -fit_w -fit_space ${space} -fit_initp ${fit_initp}
