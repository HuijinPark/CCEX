#!/bin/sh


variable_name="PULSE"
variable_values="1 4 12 64 128"

data_directory="./data/"
data_filename_format="CCE2_DiaP1_1ppm_Pulse_PULSE_Order_2_conv400"

result="./data/CCE2_DiaP1_1ppm_Pulse_PULSE_Order_2_intensity"

# Intensity plot
interpolation_x=10
interpolation_y=10
python CoherenceAnalyzer.v1/run.py -d "${data_directory}" -fi "${data_filename_format}" -vn "${variable_name}" -v ${variable_values} -uc "ns" -pl_i -interpol_x 1 -interpol_y 1 -fo "${result}" 

