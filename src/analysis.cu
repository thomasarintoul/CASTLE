__global__ void process_pairs_kernel(std::size_t pair_start, std::size_t pair_end,
                                     const std::vector<double>* input_data, // TODO: don't use host side vector
                                     std::size_t N_part,
                                     const int N_components,
                                     const int N_histograms,
                                     const std::vector<InputProperty>* properties, // TODO: don't use host side vector
                                     std::vector<double>* output_data, // TODO: don't use host side vector
                                     std::vector<OutputChannel>* output_channels) // TODO: don't use host side vector
{
    

    // Get thread index and stride length
    const std::size_t thread_index = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t stride = blockDim.x * gridDim.x;

    // Loop over pairs, calculate i and j from the pair index, and compute the structure function for each pair
    int last_reported_percent = 0;


    for (std::size_t pair_index = pair_start + thread_index; pair_index < pair_end; pair_index += stride) {
        
        // TODO: Remove this - can't be used here
        if (pair_index % 10000000 == 0){
            const std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::cout << "Time: " << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S")
                      << " - Processing pair index: " << pair_index << " of " << pair_end << "\n";
        }


        // Calculate i and j from the pair index
        //TODO: can't be used here - needs to return references or a device-compatible struct
        std::vector<std::size_t> indices = calc_i_j(pair_index, N_part);
        std::size_t i = indices[0];
        std::size_t j = indices[1];

        if (i >= N_part || j >= N_part || j <= i)
                continue;  // skip this iteration
        if (j >= N_part) j = N_part-1;


        // Get the positions of particles i and j
        // Note: by definition, x y z are the first three components of the input_data array
        std::array<double, 3> pos_i = get_position(input_data, N_part, i);
        std::array<double, 3> pos_j = get_position(input_data, N_part, j);


        // Compute distance and other relevant quantities between particle i and j
        const int storage_rows = N_histograms + 2;
        const int N_bins = output_data.size() / storage_rows;
        //TODO: Can't use vector like this on device
        const std::vector<double> bin_edges(output_data.begin(), output_data.begin() + N_bins);

        double distance = calc_separation(pos_i, pos_j);
        int bin_index = find_bin(distance, bin_edges, N_bins);
        if (bin_index < 0) {
            continue;
        }

        // Increment the relevant histogram count (counts is always the first "column" after bins in zeroth "column")
        //TODO: Make atomic
        output_data[N_bins * 1 + bin_index] += 1;


        // For each property-order-statistic combination, compute the
        // relevant quantity and update the corresponding histogram
        for (const OutputChannel& channel : output_channels) {
            // Get config values for prop-stat-order calculation
            std::size_t input_property_index = channel.input_property_index;
            Statistic statistic = channel.statistic;
            int order = channel.order;
            std::size_t output_histogram_index = channel.output_histogram_index;
            std::vector<size_t> input_components = properties[input_property_index].components;

            // Compute the relevant quantity based on the statistic and order
            double value = compute_statistic(input_data, N_part, i, j, input_components, statistic, order);
            
            // Add the computed value to the appropriate histogram in the output_data array
            //TODO: Make atomic
            output_data[N_bins * (2 + output_histogram_index) + bin_index] += value;
        }
    }

    __syncthreads();  // Ensure all threads have completed before exiting the kernel
}
