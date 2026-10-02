#!/bin/bash

BASE_DIR=./

snap=127
ANALYSIS_PARAM_FILE=./config_files/analysis_param.txt
haloes=("1")
properties=("T" "n_H" "P-KB" "z" "B")


POS_FILE=./test_inputs/input_pos.bin
SCALAR_FILE=./test_inputs/input_scalar_T.bin
IO_PARAM_FILE=./config_files/io_param.txt



build/sf_cuda -pos $BASE_DIR$POS_FILE \
            -scalar $BASE_DIR$SCALAR_FILE -an_param $ANALYSIS_PARAM_FILE \
            -io_param $IO_PARAM_FILE -snap 127 -prop T

