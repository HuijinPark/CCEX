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

result="./data/CCE2_DiaP1_1ppm_Pulse_${pulse}_Order_2"

######################################
# Ensemble average + plot
#
# Ensemble average option   : -ea
# Unit conversion           : -uc [unit] (use default unit (ms))
# Nan to num                : -ntn
python CoherenceAnalyzer.v1/run.py -d "${data_directory}" -fi "${data_filename_format}" -vn "${variable_name1}" "${variable_name2}" -v ${variable_values1} -v ${variable_values2} -uc "ms" -ea -ntn -fo "${result}" -pl #-vb

######################################
# Ensemble average only
#
# Ensemble average option   : -ea
# Unit conversion           : -uc [unit] (use default unit (ms))
# Nan to num                : -ntn
python CoherenceAnalyzer.v1/run.py -d "${data_directory}" -fi "${data_filename_format}" -vn "${variable_name1}" "${variable_name2}" -v ${variable_values1} -v ${variable_values2} -uc "ms" -ea -ntn -fo "${result}"

