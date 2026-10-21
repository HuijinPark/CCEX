"""The input layout (A64 configuration block, A65 blocks) for the input generators of run.sh.

The cases of run.sh name gyrofile, bathfile, nstate, seed, coordinate_frame_rotation, outfile, ... at the top of the
input and Defect.axis_file as one name, as they were written before 2026-10-15, and a case's overrides replace those keys
one by one. configuration_block() puts them into their blocks just before the input is written
(../input_layout.py nest).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from input_layout import nest   # noqa: E402


def configuration_block(cfg):
    return nest(cfg)
