#ifndef VSF_FUNCTIONS_H
    #define VSF_FUNCTIONS_H

    #include <iostream>     // strings
    #include <vector>       // vectors
    #include <algorithm>    // vector ops on arrays
    #include <string>        // strings
    #include <unordered_map> // hash maps

    // #include <vector>
    #include <cmath>
    #include <array>
    // #include <random>
    // #include <iterator>
    // #include <numeric>
    // #include <tuple>
    // #include <fstream>
    // #include "cnpy.h"
    // #include <Eigen/Dense> // For vector operations
    // #include </home/user/c23052501/local/include/tqdm.h> // progress bar
    #include <cuda_runtime.h>  // Core CUDA runtime types and functions



    std::vector<double> linspace(double start, double end, int num);


    std::vector<double> logspace(double start, double end, int num);


    std::vector<std::array<double, 3>> load_bin_3d(const std::string& filename); 


    std::vector<double> load_bin_1d(const std::string& filename);


    void save_to_binary(const std::string& filename, const std::vector<double>& data);


    std::unordered_map<std::string, std::string> read_parameters(const std::string& param_file);


    __device__ int calc_i(int N_cells, long long k_stride);


    __device__ int calc_ti(int N_cells, int i);


    __device__ int calc_j(int N_cells, long long k_stride, int i, int Ti);


    __device__ int find_bin(double separation, double* bins, int n_bins);


    __global__ void vsf_kernel(
        int N, // number of cells
        int n_bins, // number of bins for the structure function
        double* bins, // array of bin edges
        long long iter_start, // starting index for this kernel launch
        long long iter_end, // ending index for this kernel launch
        double* global_bin_hist, // global histogram for bin counts
        double3* gas_positions,  // point to position data in global memory
        double3* gas_velocities = nullptr, // point to velocity data in global memory
        double* gas_scalars = nullptr, // point to scalar data in global memory - nullptr if not used
        double* global_vsf = nullptr, // global histogram for sum of velocity differences in each bin
        double* global_vsf_parr = nullptr, // global histogram for sum of velocity differences in each bin
        double* global_vsf_perp = nullptr, // global histogram for sum of velocity differences in each bin
        double* global_scalar = nullptr // global histogram for sum of velocity differences in each bin
        );


    

#endif // VSF_FUNCTIONS_H