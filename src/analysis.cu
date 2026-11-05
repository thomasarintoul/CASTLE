#include "data.hpp"
#include "analysis.cpp"

#ifdef __CUDACC__
  #define HD __host__ __device__
  #define H  __host__
#else
  #define HD
  #define H
#endif



__global__ void process_pairs_kernel(std::size_t pair_start, std::size_t pair_end,
                                     const double* input_data,
                                     std::size_t N_part,
                                     const int N_components,
                                     std::size_t N_bins, 
                                     const int N_histograms,
                                     const DeviceProperty* properties, 
                                     double* output_data, 
                                     const OutputChannel* output_channels,
                                     const std::size_t output_channel_count) 
{
    

    // Get thread index and stride length
    const std::size_t thread_index = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t stride = blockDim.x * gridDim.x;

    // Loop over pairs, calculate i and j from the pair index, and compute the structure function for each pair
    // const int storage_rows = N_histograms + 2;
    const double* bin_edges = output_data;
    
    for (std::size_t pair_index = pair_start + thread_index; pair_index < pair_end; pair_index += stride) {
        
        

        // Calculate i and j from the pair index
        Pair_Index indices = calc_i_j(pair_index, N_part);
        std::size_t i = indices.i;
        std::size_t j = indices.j;

        if (i >= N_part || j >= N_part || j <= i)
                continue;  // skip this iteration


        // Get the positions of particles i and j
        // Note: by definition, x y z are the first three components of the input_data array
        std::array<double, 3> pos_i = get_position(input_data, N_part, i);
        std::array<double, 3> pos_j = get_position(input_data, N_part, j);


        // Compute distance and other relevant quantities between particle i and j
        double distance = calc_separation(pos_i, pos_j);
        int bin_index = find_bin(distance, bin_edges, N_bins);
        if (bin_index < 0) {
            continue;
        }

        // Increment the relevant histogram count (counts is always the first "column" after bins in zeroth "column")
        atomicAdd(&output_data[N_bins * 1 + bin_index], 1.0);


        // For each property-order-statistic combination, compute the
        // relevant quantity and update the corresponding histogram
        for (std::size_t c = 0; c < output_channel_count; ++c) {
            const OutputChannel& channel = output_channels[c];
            // Get config values for prop-stat-order calculation
            std::size_t input_property_index = channel.input_property_index;
            Statistic statistic = channel.statistic;
            int order = channel.order;
            std::size_t output_histogram_index = channel.output_histogram_index;
            const std::size_t* input_components = properties[input_property_index].components;

            // Compute the relevant quantity based on the statistic and order
            double value = compute_statistic(input_data, N_part, i, j, input_components, statistic, order);
            
            // Add the computed value to the appropriate histogram in the output_data array
            atomicAdd(&output_data[N_bins * (2 + output_histogram_index) + bin_index], value);
        }
    }

    __syncthreads();  // Ensure all threads have completed before exiting the kernel
}
