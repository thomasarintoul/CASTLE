// defines the loop for CPU computation of the structure function

output_array run_cpu_sf(
    const input_array& input_data,  // 2D array of input data (N_part x N_components)
    const int N_part,               // number of particles
    const int N_components,         // number of components in the input data
    const int N_histograms,         // number of histograms to compute
    output_array& output_data,      // 2D array to hold the computed histograms (N_histograms x N_bins)
    mapping& output_channels,       // mapping of output channels to input and output histogram indices
) {

    // From provided output_data array, get the number of bins in the histograms and the min and max separation values for the histograms
    /// Bin lower edges are in the first "column" of the output_data array, and the counts/statistics are in the subsequent "columns"
    /// The array is actually 1D contiguous, but we treat it as 2D for clarity. First N_bins values are bin edges, next is counts, etc.
    const int N_bins = output_data.values.size() / N_histograms;

    // Loop over all pairs of particles to compute structure functions
    for (int i = 0; i < N_part; ++i) {
        for (int j = i + 1; j < N_part; ++j) {
            // Get the positions of particles i and j
            /// Note that by definition, x y z are the first three components of the input_data array
            std::array<double, 3> pos_i = {input_data[N_part*0 + i], input_data[N_part*1 + i], input_data[N_part*2 + i]};
            std::array<double, 3> pos_j = {input_data[N_part*0 + j], input_data[N_part*1 + j], input_data[N_part*2 + j]};

            // Compute distance and other relevant quantities between particle i and j
            double distance = calc_separation(pos_i, pos_j);
            int bin_index = find_bin(distance, output_data);

            // Increment the relevant histogram count (counts is always the first "column" after bins in zeroth "column")
            output_data[N_histograms + bin_index] += 1;


            // For each property-order-statistic combination, compute the relevant quantity and update the corresponding histogram
            
            
        }
    }

    return output;
}