#!/bin/sh

# File names
# variable_name will be replaced with the variable_values
variable_name1="CONF"
variable_values1="$(seq -s ' ' 1 1 20)"

data_directory="./results/"
data_filename_format="CCE2_Diamond_NV_natab_CONF_wiDiv"
data_directory_ref="./results/"
data_filename_format_ref="CCE2_Diamond_NV_natab_CONF_noDiv"

result="./results/CCE2_Diamond_NV_natab_CONF_ErrCorr"

######################################
# Error correction and save the corrected data
CodePath="$(git rev-parse --show-toplevel)/example/z.Analyzer/CoherenceAnalyzer.v1/run.py"
python ${CodePath} -d "${data_directory}" -fi "${data_filename_format}" -err_d "${data_directory_ref}" -err_fi "${data_filename_format_ref}" -vn "${variable_name1}" -v ${variable_values1}  -uc "ms" -err -ntn -fo "${result}" #-vb


