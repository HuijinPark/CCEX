#!/bin/sh

# File names
# variable_name will be replaced with the variable_values
variable_name1="PREFIX"
variable_values1="single multi2"
#variable_values1="hf3_qd0 hf3_qd2"

data_directory="./results/"
data_filename_format="CCE2_Diamond_NV_natab_PREFIX_1_wiDiv"

result="./CCE2_Diamond_NV_natab_Comp_Multi"
#result="./CCE2_Diamond_NV_natab_Comp_Qd"

CodePath="$(git rev-parse --show-toplevel)/example/z.Analyzer/CoherenceAnalyzer.v1/run.py"
######################################
# Ensemble average + plot
#
# Ensemble average option   : -ea
# Unit conversion           : -uc [unit] (use default unit (ms))
# Nan to num                : -ntn
python ${CodePath} -d "${data_directory}" -fi "${data_filename_format}" -vn "${variable_name1}" -v ${variable_values1} -uc "us" -ntn -fo "${result}" -pl_m -ylim 0 1 #-vb

