#!/bin/sh

variable_name="PULSE"
variable_values="1 4 12 64 128"

data_directory="./data/"
data_filename_format="CCE2_DiaP1_1ppm_Pulse_PULSE_Order_2_conv400"

result="./data/CCE2_DiaP1_1ppm_Pulse_PULSE_Order_2"

###############################
# Single plot
python CoherenceAnalyzer.v1/run.py -d "${data_directory}" -fi "${data_filename_format}" -vn "${variable_name}" -v ${variable_values} -uc "ns" -pl -ntn -fo "${result}" -ylim 0 1

###############################
# Single plot with ftting
#
# Fitting : -fit
result="./data/CCE2_DiaP1_1ppm_Pulse_PULSE_Order_2_Fitting"
python CoherenceAnalyzer.v1/run.py -d "${data_directory}" -fi "${data_filename_format}" -vn "${variable_name}" -v ${variable_values} -uc "ns" -pl -ntn -fo "${result}" -ylim 0 1 -fit
