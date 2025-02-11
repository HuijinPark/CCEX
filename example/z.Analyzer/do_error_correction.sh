#!/bin/sh

# File names
# variable_name will be replaced with the variable_values
variable_name1="CONF"
variable_name2="STATE"
variable_values1="$(seq -s ' ' 1 1 20)"
variable_values2="$(seq -s ' ' 1 1 20)"

pulse="128"
data_directory="./data/CCE2_1ppm_Pulse_${pulse}/"
data_filename_format="CCE2_DiaP1_1ppm_Pulse_${pulse}_Order_2_wD_CONF{STATE}"
data_directory_ref="./data/CCE2_1ppm_Pulse_${pulse}/"
data_filename_format_ref="CCE2_DiaP1_1ppm_Pulse_${pulse}_Order_2_nD_CONF{STATE}"

result="./data/CCE2_1ppm_Pulse_${pulse}/CCE2_DiaP1_1ppm_Pulse_${pulse}_Order_2_EC_CONF{STATE}"

######################################
# Error correction and save the corrected data
python CoherenceAnalyzer.v1/run.py -d "${data_directory}" -fi "${data_filename_format}" -err_d "${data_directory_ref}" -err_fi "${data_filename_format_ref}" -vn "${variable_name1}" "${variable_name2}" -v ${variable_values1} -v ${variable_values2} -uc "ms" -err -ntn -fo "${result}" #-vb


