#!/bin/sh

# File names
# variable_name will be replaced with the variable_values
variable_name1="ZFSIO"
variable_values1="with without"

data_directory="./0G_ZFSIO_ZFS/rawdata/"
data_filename_format="gCCE2_DiaP1_10ppm_B0_0G_1_wiDiv"

result="./Coherence_Comp_wi_wo_ZFS"

CodePath="$(git rev-parse --show-toplevel)/example/z.Analyzer/CoherenceAnalyzer.v1/run.py"
######################################
# Ensemble average + plot
#
# Ensemble average option   : -ea
# Unit conversion           : -uc [unit] (use default unit (ms))
# Nan to num                : -ntn
python ${CodePath} -d "${data_directory}" -fi "${data_filename_format}" -vn "${variable_name1}" -v ${variable_values1} -uc "us" -ntn -fo "${result}" -pl_m #-vb

