#!/bin/sh

variable_name="PULSE"
variable_values="1 4 12 64 128"

data_directory="./data/"
data_filename_format="CCE2_DiaP1_1ppm_Pulse_PULSE_Order_2_conv400"

result="./data/CCE2_DiaP1_1ppm_Pulse_PULSE_Order_2"

###############################
# MultiPlot
python CoherenceAnalyzer.v1/run.py -d "${data_directory}" -fi "${data_filename_format}" -vn "${variable_name}" -v ${variable_values} -uc "ns" -pl_m -ntn -fo "${result}" -ylim 0 1

