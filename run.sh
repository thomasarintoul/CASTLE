#!/bin/bash

scalar_prop="T" 
BASE_DIR=/home/user/c23052501/CASTLE/TSF_subset/
ANALYSIS_PARAM_FILE=$BASE_DIR/cuda_vsf_param.txt

echo ''
echo ''
echo '// STARTING VSF CALCULATION //'
echo '- Model: TNG'

./sf_cuda.out -pos $BASE_DIR/sur5_2kpc_tng_subset_pos.bin \
               -scalar $BASE_DIR/sur5_2kpc_tng_subset_T.bin -an_param $ANALYSIS_PARAM_FILE \
               -io_param $BASE_DIR/io_tng.txt -snap 127 -prop $scalar_prop

echo ''
echo ''
echo '// STARTING VSF CALCULATION //'
echo '- Model: WTC'

./sf_cuda.out -pos $BASE_DIR/sur5_2kpc_wtc_subset_pos.bin \
               -scalar $BASE_DIR/sur5_2kpc_wtc_subset_T.bin -an_param $ANALYSIS_PARAM_FILE \
               -io_param $BASE_DIR/io_wtc.txt -snap 127 -prop $scalar_prop


echo ''
echo ''
echo '// STARTING VSF CALCULATION //'
echo '- Model: STC'

./sf_cuda.out -pos $BASE_DIR/sur5_2kpc_stc_subset_pos.bin \
               -scalar $BASE_DIR/sur5_2kpc_stc_subset_T.bin -an_param $ANALYSIS_PARAM_FILE \
               -io_param $BASE_DIR/io_stc.txt -snap 127 -prop $scalar_prop


echo ''
echo ''
echo '// STARTING VSF CALCULATION //'
echo '- Model: TNG - Hot'

./sf_cuda.out -pos $BASE_DIR/sur5_2kpc_tng_hot_subset_pos.bin \
               -scalar $BASE_DIR/sur5_2kpc_tng_hot_subset_T.bin -an_param $ANALYSIS_PARAM_FILE \
               -io_param $BASE_DIR/io_tng_hot.txt -snap 127 -prop $scalar_prop


echo ''
echo ''
echo '// STARTING VSF CALCULATION //'
echo '- Model: WTC - Hot'

./sf_cuda.out -pos $BASE_DIR/sur5_2kpc_wtc_hot_subset_pos.bin \
               -scalar $BASE_DIR/sur5_2kpc_wtc_hot_subset_T.bin -an_param $ANALYSIS_PARAM_FILE \
               -io_param $BASE_DIR/io_wtc_hot.txt -snap 127 -prop $scalar_prop


echo ''
echo ''
echo '// STARTING VSF CALCULATION //'
echo '- Model: STC - Hot'

./sf_cuda.out -pos $BASE_DIR/sur5_2kpc_stc_hot_subset_pos.bin \
               -scalar $BASE_DIR/sur5_2kpc_stc_hot_subset_T.bin -an_param $ANALYSIS_PARAM_FILE \
               -io_param $BASE_DIR/io_stc_hot.txt -snap 127 -prop $scalar_prop
