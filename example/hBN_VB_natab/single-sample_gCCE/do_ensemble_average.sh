#!/bin/sh

# File names
# variable_name will be replaced with the variable_values
variable_name1="CONF"
variable_values1="$(seq -s ' ' 1 1 20)"
variable_name2="STAT"
variable_values2="$(seq -s ' ' 1 1 20)"

data_directory="./results/"
data_filename_format="CCE2_Diamond_NV_natab_CONF_stateSTAT_ErrCorr"

result="./results/CCE2_Diamond_NV_natab_EnsAvg"

CodePath="$(git rev-parse --show-toplevel)/example/z.Analyzer/CoherenceAnalyzer.v1/run.py"
######################################
# Ensemble average + plot
#
# Ensemble average option   : -ea
# Unit conversion           : -uc [unit] (use default unit (ms))
# Nan to num                : -ntn
python ${CodePath} -d "${data_directory}" -fi "${data_filename_format}" -vn "${variable_name1}" "${variable_name2}"  -v ${variable_values1} -v ${variable_values2} -uc "ms" -ea -ntn -fo "${result}" -pl -fit -fit_initp 2 -fit_space 5 #-vb

