#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>
#include <algorithm>
#include <cstdint>
#include "toml.hpp"
#include "analysis.cu"



// Defines the host-side orchestration for CUDA computation of the structure function.
// This function validates inputs, prepares device buffers, launches the kernel, and
// combines the per-GPU results.
void run_cuda_sf(
    const std::vector<double>& input_data,              // 2D array of input data (N_part x N_components)
    std::size_t N_part,                           // number of particles
    const int N_components,                     // number of components in the input data
    const int N_histograms,                     // number of histograms to compute
    const std::vector<InputProperty>& properties,                  // vector of InputProperty structs defining the properties to compute
    std::vector<double>& output_data,                  // 2D array to hold the computed histograms (N_histograms x N_bins)
    std::vector<OutputChannel>& output_channels,  // mapping of output channels to input and output histogram indices
    int N_gpus = 1                                // number of GPUs to use (default is 1)  
) {

    if (N_histograms < 0 || output_data.empty()) {
        throw std::runtime_error("Invalid histogram configuration.");
    }

    if (N_components < 3) {
        throw std::runtime_error("At least three input components are required for positions.");
    }
    const std::size_t required_input_values =
        N_part * static_cast<std::size_t>(N_components);
    if (input_data.size() < required_input_values) {
        throw std::runtime_error("Input data is smaller than the configured dimensions.");
    }

    // From provided output_data array, get the number of bins in the histograms and the min and max separation values for the histograms
    /// Bin lower edges are in the first "column" of the output_data array, and the counts/statistics are in the subsequent "columns"
    /// The array is actually 1D contiguous, but we treat it as 2D for clarity. First N_bins values are bin edges, next is counts, etc.
    const int storage_rows = N_histograms + 2;
    if (storage_rows <= 0 || output_data.size() % static_cast<std::size_t>(storage_rows) != 0) {
        throw std::runtime_error("Output data size does not match the histogram configuration.");
    }
    const int N_bins = output_data.size() / storage_rows;
    if (N_bins < 2) {
        throw std::runtime_error("At least two bin edges are required.");
    }

    for (const OutputChannel& channel : output_channels) {
        if (channel.input_property_index >= properties.size() ||
            channel.output_histogram_index >= static_cast<std::size_t>(N_histograms)) {
            throw std::runtime_error("Output channel contains an invalid index.");
        }
        const std::size_t required_components =
            channel.statistic == Statistic::ScalarDifference ? 1 : 3;
        if (properties[channel.input_property_index].components.size() < required_components) {
            throw std::runtime_error("Statistic requires more input components than configured.");
        }
    }
    const std::vector<double> bin_edges(output_data.begin(), output_data.begin() + N_bins);
    std::cout << "Number of bins: " << N_bins << "\n";
    std::cout << "Bin edges: ";
    for (const auto& edge : bin_edges) {
        std::cout << edge << " ";
    }
    std::cout << "\n";

    // Convert InputProperty to DeviceProperty for device-side use
    std::vector<DeviceProperty> device_properties = conv_device_prop(properties);
    
    // TODO: Remove this testing limit later.
    N_part = std::min<std::size_t>(N_part, 100000);
    
    
    // Total number of unique pairs (i, j) with i < j: N * (N - 1) / 2.
    // Use 64-bit integers because this overflows int quickly (N = 100,000 gives ~5 billion).
    const std::uint64_t num_particles = static_cast<std::uint64_t>(N_part);
    const std::uint64_t total_pairs   = num_particles * (num_particles - 1) / 2;
    // If N_gpus > 1, we will divide the total pairs into chunks for each GPU to process
    std::vector<std::vector<double>> per_gpu_output_data;
    if (N_gpus > 1) {
        // Create local output arrays for each GPU to avoid race conditions
        // Each local array is initialised to zero
        per_gpu_output_data.assign(N_gpus, std::vector<double>(output_data.size(), 0.0));
        
        // Copy bin edges from the main output_data array to each GPU's local output array
        for (int gpu_id = 0; gpu_id < N_gpus; ++gpu_id) {
            std::copy(output_data.begin(), output_data.begin() + N_bins, per_gpu_output_data[gpu_id].begin());
        }
    } else {
        // If only one GPU, we can use the main output_data array directly
        per_gpu_output_data.assign(1, output_data);
    }



    #pragma omp parallel num_threads(N_gpus)
    {
        const int gpu_id = omp_get_thread_num();
        cudaSetDevice(gpu_id);  // Set the current thread to use the correct GPU
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "CUDA error on GPU " << gpu_id << ": " << cudaGetErrorString(err) << std::endl;
        }


        double* d_input_data = nullptr;
        double* d_output_data = nullptr;
        // double* d_properties = nullptr;
        DeviceProperty* d_device_properties = nullptr;
        OutputChannel* d_output_channels = nullptr;

        cudaMalloc(&d_input_data, input_data.size() * sizeof(double));
        cudaMalloc(&d_output_data, per_gpu_output_data[gpu_id].size() * sizeof(double));
        cudaMalloc(&d_device_properties, device_properties.size() * sizeof(DeviceProperty));
        cudaMalloc(&d_output_channels, output_channels.size() * sizeof(OutputChannel));

        cudaMemcpy(d_input_data, 
                   input_data.data(), 
                   input_data.size() * sizeof(double), 
                   cudaMemcpyHostToDevice);
        cudaMemcpy(d_output_data, 
                   per_gpu_output_data[gpu_id].data(), 
                   per_gpu_output_data[gpu_id].size() * sizeof(double), 
                   cudaMemcpyHostToDevice);
        cudaMemcpy(d_device_properties, 
                   device_properties.data(), 
                   device_properties.size() * sizeof(DeviceProperty), 
                   cudaMemcpyHostToDevice);
        cudaMemcpy(d_output_channels, 
                   output_channels.data(), 
                   output_channels.size() * sizeof(OutputChannel), 
                   cudaMemcpyHostToDevice);
        std::cout << "GPU " << gpu_id << ": copied data to device memory.\n";


        // Kernel Config (per GPU)

        int blockSize = 256;  // Number of threads per block
        int numSM;
        cudaDeviceGetAttribute(&numSM, cudaDevAttrMultiProcessorCount, gpu_id);
        // TODO: Benchmark with different grid size factors
        int gridSize = numSM * 8;  // Number of blocks in the grid

        // Calculate number of pairs per GPU (excl. remainder)
        const std::size_t pairs_per_gpu = total_pairs / N_gpus;
        // Calculate the remainder pairs that need to be distributed among the first few GPUs
        const std::size_t remainder_pairs = total_pairs % N_gpus;
        
        const std::size_t pair_start = gpu_id * pairs_per_gpu + \
                                       std::min<std::size_t>(gpu_id, remainder_pairs);
        const std::size_t pair_end = pair_start + pairs_per_gpu + \
                                     (gpu_id < remainder_pairs ? 1 : 0);


        std::cout << "GPU " << gpu_id << " processing pairs from " << pair_start << " to " << pair_end << "\n";


        // Launch the kernel to process pairs on the GPU
        process_pairs_kernel<<<gridSize, blockSize>>>(
            pair_start,
            pair_end,
            d_input_data,
            N_part,
            N_components,
            N_bins, 
            N_histograms,
            d_device_properties,
            d_output_data,
            d_output_channels,
            output_channels.size()
        );

        // Check for kernel launch errors
        err = cudaGetLastError();
        if (err != cudaSuccess)
            std::cerr << "GPU " << gpu_id << " kernel launch error: " << cudaGetErrorString(err) << "\n";


        cudaDeviceSynchronize();

        err = cudaGetLastError();
        if (err != cudaSuccess)
            std::cerr << "GPU " << gpu_id << " kernel execution error: " << cudaGetErrorString(err) << "\n";

        // Copy back the per-GPU output array 
        cudaMemcpy(per_gpu_output_data[gpu_id].data(),
                   d_output_data,
                   per_gpu_output_data[gpu_id].size() * sizeof(double),
                   cudaMemcpyDeviceToHost);

        cudaFree(d_input_data);
        cudaFree(d_output_data);
        cudaFree(d_device_properties);
        cudaFree(d_output_channels);


    }



    // After processing, combine results of all threads into the main output_data array
    // In the case of single threaded, this effectively copies the results from thread_output_data[0] to output_data
    // For multi-threaded, this sums the results from all threads into the main output_data array
    for (int gpu_id = 0; gpu_id < N_gpus; ++gpu_id) {
        // Loop over all rows except the bin edges (row 0), which are already handled
        for (std::size_t idx = N_bins; idx < output_data.size(); ++idx) {
            output_data[idx] += per_gpu_output_data[gpu_id][idx];
        }
    }

}
