#include <iostream>
#include <vector>
#include <random>
#include <iterator>
#include <numeric>
#include <tuple>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <atomic>
#include <array>
#include <thread>
#include <sstream>
#include <string>
#include <unordered_map>
#include <cuda_runtime.h>  // Core CUDA runtime types and functions
#include <device_launch_parameters.h> // Needed for blockIdx, threadIdx, etc.
#include <cmath> // For sqrt and pow
#include <cstddef> // For atomic operations
#include <cstdint> // For fixed-width integer types
#include <algorithm> // min/max for clamping indices
#include <math.h> // For sqrt and pow


std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> linspaced;
    double delta = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        linspaced.push_back(start + delta * i);
    }
    return linspaced;
}


std::vector<double> logspace(double start, double end, int num) {
    std::vector<double> logspaced;
    double delta = (log10(end) - log10(start)) / (num - 1);
    for (int i = 0; i < num; ++i) {
        logspaced.push_back(pow(10, log10(start) + delta * i));
    }
    return logspaced;
}


std::vector<std::array<double, 3>> load_bin_3d(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    size_t num_bytes = file.tellg();
    file.seekg(0, std::ios::beg);

    if (num_bytes % (sizeof(double) * 3) != 0) {
        throw std::runtime_error("File size is not a multiple of 3 * sizeof(double)");
    }

    size_t N = num_bytes / (sizeof(double) * 3);
    std::vector<std::array<double, 3>> data(N);

    file.read(reinterpret_cast<char*>(data.data()), num_bytes);
    return data;
}


std::vector<double> load_bin_1d(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    size_t num_bytes = file.tellg();
    file.seekg(0, std::ios::beg);

    if (num_bytes % (sizeof(double) * 1) != 0) {
        throw std::runtime_error("File size is not a multiple of 1 * sizeof(double)");
    }

    size_t N = num_bytes / (sizeof(double) * 1);
    std::vector<double> data(N);

    file.read(reinterpret_cast<char*>(data.data()), num_bytes);
    return data;
}



void save_to_binary(const std::string& filename, const std::vector<double>& data) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }
    
    file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(double));
    file.close();
    std::cout << "Saved data to " << filename << std::endl;
}


std::unordered_map<std::string, std::string> read_parameters(const std::string& param_file) {
    // Create the unordered_map to store parameters (similar to a dict in python)
    std::unordered_map<std::string, std::string> params;
    
    // Open the parameter file
    std::ifstream file(param_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open parameter file: " << param_file << std::endl;
        std::exit(1);
    }

    // Read the file line by line
    std::string line;
    while (std::getline(file, line)) {
        // Skip lines starting with '%'
        if (!line.empty() && line[0] == '%') {
            continue;
        }

        std::istringstream iss(line);
        std::string key, value;
        if (iss >> key >> value) {
            params[key] = value;
        }
    }
    file.close();
    return params;
}


// __device__ int calc_i(int N_cells, long long k_stride){
//     // for given N_cells and k_stride, calculate i using 
//     // the formula derived from the triangular number sequence
//     return static_cast<int>(((2.0*N_cells - 1.0) - sqrt((2.0*N_cells - 1.0)*(2.0*N_cells - 1.0) - 8.0*k_stride)) / 2);
// }

__device__ int calc_i(int N_cells, long long k) {
    double a = 2.0 * (double)N_cells - 1.0;
    double val = a*a - 8.0 * (double)k;
    int i = (int)((a - sqrt(val)) * 0.5);
    return i;
}

__device__ int calc_ti(int N_cells, int i){
    // for given N_cells and row i, calculate Ti using 
    // the formula derived from the triangular number sequence
    return static_cast<int>(i*N_cells - (i*(i+1)) / 2);
}

__device__ int calc_j(int N_cells, long long k_stride, int i, int Ti){
    // for given N_cells, k_stride, i, and Ti, calculate row entry j using 
    // the formula derived from the triangular number sequence
    return static_cast<int>(i + 1 + k_stride - Ti);
}

__device__ int find_bin(double separation, double* bins, int n_bins) {

    if (separation < bins[0] || separation >= bins[n_bins - 1]) return -1;

    int left = 0;
    int right = n_bins - 1;

    while (right - left > 1) {
        int mid = (left + right) / 2;

        if (separation >= bins[mid])
            left = mid;
        else
            right = mid;
    }

    return left;
}


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
    ) {
        // get shared memory data
        extern __shared__ double shared[];

        // local histograms
        double* local_hist      = shared;
        double* local_vsf       = (global_vsf)       ? &shared[n_bins * 1] : nullptr; // allocate local_vsf only if global_vsf is provided
        double* local_vsf_parr  = (global_vsf_parr)  ? &shared[n_bins * 2] : nullptr; // * 2 to avoid overlap with local_vsf
        double* local_vsf_perp  = (global_vsf_perp)  ? &shared[n_bins * 3] : nullptr; // * 3 to avoid overlap with local_vsf and local_vsf_parr
        double* local_scalar    = (global_scalar)    ? &shared[n_bins * 4] : nullptr; // * 4 to avoid overlap with local_vsf, local_vsf_parr, and local_vsf_perp

        if (gas_velocities == nullptr && gas_scalars == nullptr) {
            printf("Error: No velocity or scalar data provided to kernel.\n");
            return;
        }
        

        // initialise local histograms to zero
        for (int i = threadIdx.x; i < n_bins; i += blockDim.x) {
            local_hist[i] = 0.0;
            if (local_vsf)      local_vsf[i] = 0.0; 
            if (local_vsf_parr) local_vsf_parr[i] = 0.0;
            if (local_vsf_perp) local_vsf_perp[i] = 0.0;
            if (local_scalar)   local_scalar[i] = 0.0;
        }

        __syncthreads(); // ensure all threads have initialised their local histograms

        // get global thread ID - thread index in block + block index * block size
        long long global_id = threadIdx.x + blockIdx.x * blockDim.x;

        // set stride length for parallel processing of pairs
        long long stride_length = blockDim.x * gridDim.x; // total number of threads in the grid

        for (long long k = iter_start + global_id; k < iter_end; k += stride_length) {
            int i = calc_i(N, k);
            int Ti = calc_ti(N, i);
            int j = calc_j(N, k, i, Ti);

            // Safety check
            if (i < 0 || i >= N || j < 0 || j >= N || j <= i)
                continue;  // skip this iteration

            if (j >= N) j = N-1;

            // load positions and velocities for the pair (i, j)
            double3 pos_i = gas_positions[i];
            double3 pos_j = gas_positions[j];

            double3 vel_i = (gas_velocities != nullptr) ? gas_velocities[i] : make_double3(0.0, 0.0, 0.0);
            double3 vel_j = (gas_velocities != nullptr) ? gas_velocities[j] : make_double3(0.0, 0.0, 0.0);

            double scalar_i = (gas_scalars != nullptr) ? gas_scalars[i] : 0.0;
            double scalar_j = (gas_scalars != nullptr) ? gas_scalars[j] : 0.0;


            // calculate separation and velocity difference for the pair
            double x_sep = pos_i.x - pos_j.x;
            double y_sep = pos_i.y - pos_j.y;
            double z_sep = pos_i.z - pos_j.z;
            double separation = sqrt(x_sep * x_sep + y_sep * y_sep + z_sep * z_sep);

            // Ignore pairs outside the configured separation range.
            if (separation < bins[0] || separation >= bins[n_bins - 1])
                continue;

            // Only calculate velocity differences if velocity data is provided
            double x_vdiff = (gas_velocities != nullptr) ? vel_i.x - vel_j.x : 0.0;
            double y_vdiff = (gas_velocities != nullptr) ? vel_i.y - vel_j.y : 0.0;
            double z_vdiff = (gas_velocities != nullptr) ? vel_i.z - vel_j.z : 0.0;
            double v_diff = (gas_velocities != nullptr) ? sqrt(x_vdiff * x_vdiff + y_vdiff * y_vdiff + z_vdiff * z_vdiff) : 0.0;


            // Only calculate scalar differences if scalar data is provided
            double scalar_diff = (gas_scalars != nullptr) ? ((scalar_i - scalar_j) * (scalar_i - scalar_j)) : 0.0; // assuming scalar is stored as a single value

            // If parallel and perpendicular projected velocity differences are requested,
            // calculate them using the separation vector and determined normal vector
            bool do_rhat = (global_vsf_parr != nullptr || global_vsf_perp != nullptr);
            bool do_nhat = (global_vsf_perp != nullptr);

            double r_hat_x = (do_rhat) ? x_sep / separation : 0.0;
            double r_hat_y = (do_rhat) ? y_sep / separation : 0.0;
            double r_hat_z = (do_rhat) ? z_sep / separation : 0.0;

            double n_hat_x = (do_nhat) ? -r_hat_y : 0.0; // arbitrary choice of normal vector in 2D plane
            double n_hat_y = (do_nhat) ? r_hat_x : 0.0;
            double n_hat_z = (do_nhat) ? 0.0 : 0.0;

            double v_parr_diff = (global_vsf_parr != nullptr) ? sqrt((x_vdiff * r_hat_x) * (x_vdiff * r_hat_x) + 
                                                              (y_vdiff * r_hat_y) * (y_vdiff * r_hat_y) + 
                                                              (z_vdiff * r_hat_z) * (z_vdiff * r_hat_z)) : 0.0;
            double v_perp_diff = (global_vsf_perp != nullptr) ? sqrt((x_vdiff * n_hat_x) * (x_vdiff * n_hat_x) + 
                                                              (y_vdiff * n_hat_y) * (y_vdiff * n_hat_y) + 
                                                              (z_vdiff * n_hat_z) * (z_vdiff * n_hat_z)) : 0.0;



            // add to local histograms
            int bin_index = find_bin(separation, bins, n_bins);

            if (bin_index >= 0 && bin_index < n_bins) {
                atomicAdd(&local_hist[bin_index], 1.0); // increment count for this bin
                if (global_vsf != nullptr) { atomicAdd(&local_vsf[bin_index], v_diff); }
                if (global_vsf_parr != nullptr) { atomicAdd(&local_vsf_parr[bin_index], v_parr_diff); }
                if (global_vsf_perp != nullptr) { atomicAdd(&local_vsf_perp[bin_index], v_perp_diff); }
                if (global_scalar != nullptr) { atomicAdd(&local_scalar[bin_index], scalar_diff); }
            }
        }

        // wait for all threads to finish processing pairs
        __syncthreads();

        // aggregate local histograms into global histograms
        for (int i = threadIdx.x; i < n_bins; i += blockDim.x){
            atomicAdd(&global_bin_hist[i], local_hist[i]);
            if (global_vsf != nullptr) { atomicAdd(&global_vsf[i], local_vsf[i]); }
            if (global_vsf_parr != nullptr) { atomicAdd(&global_vsf_parr[i], local_vsf_parr[i]); }
            if (global_vsf_perp != nullptr) { atomicAdd(&global_vsf_perp[i], local_vsf_perp[i]); }
            if (global_scalar != nullptr) { atomicAdd(&global_scalar[i], local_scalar[i]); }
        }
    }

