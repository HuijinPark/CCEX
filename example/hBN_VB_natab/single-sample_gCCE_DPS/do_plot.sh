#!/bin/sh

# File names
# variable_name will be replaced with the variable_values
variable_name1="PREFIX"
variable_values1="wi_OH wo_OH"
#variable_values1="hf3_qd0 hf3_qd2"

data_directory="./results/"
data_filename_format="CCE2_hBN_VB_PREFIX_1_state1_wiDiv"

result="./CCE2_hBN_VB_wi_wo_OH"
#result="./CCE2_hBN_VB_Comp_Qd"

CodePath="$(git rev-parse --show-toplevel)/example/z.Analyzer/CoherenceAnalyzer.v1/run.py"
######################################
# Ensemble average + plot
#
# Ensemble average option   : -ea
# Unit conversion           : -uc [unit] (use default unit (ms))
# Nan to num                : -ntn
python ${CodePath} -d "${data_directory}" -fi "${data_filename_format}" -vn "${variable_name1}" -v ${variable_values1} -uc "ns" -ntn -fo "${result}" -pl_m #-vb

