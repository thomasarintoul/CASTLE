#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
#include <stdexcept>
#include <string>
#include <vector>
#include "toml.hpp"
#include "io.cpp"
#include "data.cpp"
#include "config.cpp"
#include "analysis.cpp"


// defines the loop for CPU computation of the structure function

void run_cpu_sf(
    const std::vector<double>& input_data,              // 2D array of input data (N_part x N_components)
    int N_part,                           // number of particles
    const int N_components,                     // number of components in the input data
    const int N_histograms,                     // number of histograms to compute
    std::vector<InputProperty> properties,                  // vector of InputProperty structs defining the properties to compute
    std::vector<double>& output_data,                  // 2D array to hold the computed histograms (N_histograms x N_bins)
    std::vector<OutputChannel> output_channels  // mapping of output channels to input and output histogram indices
) {

    // From provided output_data array, get the number of bins in the histograms and the min and max separation values for the histograms
    /// Bin lower edges are in the first "column" of the output_data array, and the counts/statistics are in the subsequent "columns"
    /// The array is actually 1D contiguous, but we treat it as 2D for clarity. First N_bins values are bin edges, next is counts, etc.
    const int storage_rows = N_histograms + 2;
    const int N_bins = output_data.size() / storage_rows;
    const std::vector<double> bin_edges(output_data.begin(), output_data.begin() + N_bins);
    std::cout << "Number of bins: " << N_bins << "\n";
    std::cout << "Bin edges: ";
    for (const auto& edge : bin_edges) {
        std::cout << edge << " ";
    }
    std::cout << "\n";


    // Total number of unique pairs (i, j) with i < j: N * (N - 1) / 2.
    // Use 64-bit integers because this overflows int quickly (N = 100,000 gives ~5 billion).
    const std::uint64_t num_particles = static_cast<std::uint64_t>(N_part);
    const std::uint64_t total_pairs   = num_particles * (num_particles - 1) / 2;
    int last_reported_percent = -1;  // -1 so the first iteration always prints


    //TODO: Remove this later
    N_part = 1000; // for testing, limit to 1000 particles to avoid long runtimes


    // Loop over all pairs of particles to compute structure functions
    for (int i = 0; i < N_part; ++i) {


       // Pairs already finished before starting row i.
        // Row k contains (N - 1 - k) pairs, so summing rows 0..i-1 gives
        // i*N - i*(i+1)/2.
        const std::uint64_t row_index    = static_cast<std::uint64_t>(i);
        const std::uint64_t pairs_done   = row_index * num_particles
                                        - row_index * (row_index + 1) / 2;
        
        const int percent_done = total_pairs
            ? static_cast<int>(100.0 * pairs_done / total_pairs)
            : 100;  // guards against divide-by-zero when N_part < 2


        // Only print when the whole-number percentage changes, to avoid flooding the output
        if (percent_done != last_reported_percent) {
            last_reported_percent = percent_done;
            const auto current_time = std::chrono::system_clock::to_time_t(
                std::chrono::system_clock::now());
            std::cout << "Processed " << pairs_done << " of " << total_pairs
                    << " pairs (" << percent_done << "%) at "
                    << std::ctime(&current_time);  // ctime appends its own newline
        }


        for (int j = i + 1; j < N_part; ++j) {
            // Get the positions of particles i and j
            /// Note that by definition, x y z are the first three components of the input_data array
            std::array<double, 3> pos_i = {input_data[N_part*0 + i], input_data[N_part*1 + i], input_data[N_part*2 + i]};
            std::array<double, 3> pos_j = {input_data[N_part*0 + j], input_data[N_part*1 + j], input_data[N_part*2 + j]};


            // Compute distance and other relevant quantities between particle i and j
            double distance = calc_separation(pos_i, pos_j);
            int bin_index = find_bin(distance, bin_edges, N_bins);
            if (bin_index < 0) {
                continue;
            }


            // Increment the relevant histogram count (counts is always the first "column" after bins in zeroth "column")
            output_data[N_bins * 1 + bin_index] += 1;


            // For each property-order-statistic combination, compute the 
            // relevant quantity and update the corresponding histogram
            for (const OutputChannel& channel : output_channels) {
                // get the input property index for this channel
                std::size_t input_property_index = channel.input_property_index;

                // get statistic and order for this channel
                Statistic statistic = channel.statistic;
                int order = channel.order;

                // get the output histogram index for this channel
                std::size_t output_histogram_index = channel.output_histogram_index;

                // get the input components for this channel
                std::vector<size_t> input_components = properties[input_property_index].components;


                // Compute the relevant quantity based on the statistic and order
                double value = compute_statistic(input_data, N_part, i, j, input_components, statistic, order);
                
                // Add the computed value to the appropriate histogram in the output_data array
                output_data[N_bins * (2 + output_histogram_index) + bin_index] += value;
            }
            
            
        }
    }
}

int main(int argc, char* argv[]) {    
    // Call the function to configure the analysis and retrieve the necessary data structures
    Core_Data core_data;
    try{
        core_data = configure_analysis(argc, argv);
    }
    catch (const std::runtime_error& e) {
        std::cerr << "Configuration error: " << e.what() << "\n";
        return 1;
    }


    // Call the function to run the CPU structure function computation
    run_cpu_sf(core_data.input_data,         // 2D array of input data (N_part x N_components)
               core_data.N_part,             // N of elements in the input data (number of particles)
               core_data.N_components,       // N of input data components (3 for vectors, 1 for scalars)
               core_data.N_histograms,       // N of output hists (i.e. number stat-order-prop combinations)
               core_data.properties,         // vector of InputProperty structs defining input properties
               core_data.output_data,        // 2D array to hold the computed histograms (N_histograms x N_bins)
               core_data.output_channels);   // mapping of output channels to input and output histogram indices


    // Print the first 10 values of each row of the output_data array for verification
    std::cout << "\n========== Output Data Summary ==========" << "\n";
    std::cout << "Output data (first 10 values of each histogram):\n";

    // for each histogram, print the first 10 values of the corresponding row in the output_data array
    
    // start with bin edges and counts, which are always the first two rows in the output_data array
    std::cout << ">> Bin edges (row 0): ";
    for (int bin = 0; bin < std::min(10, static_cast<int>(core_data.output_data.size() / (core_data.N_histograms + 2))); ++bin) {
        std::cout << core_data.output_data[bin] << " ";
    }
    std::cout << "\n";
    std::cout << ">> Counts (row 1): ";
    for (int bin = 0; bin < std::min(10, static_cast<int>(core_data.output_data.size() / (core_data.N_histograms + 2))); ++bin) {
        std::cout << core_data.output_data[bin + (1 * (core_data.output_data.size() / (core_data.N_histograms + 2)))] << " ";
    }
    std::cout << "\n";

    // then loop over the property-statistic-order combinations, which are stored in the subsequent rows of the output_data array
    for (int hist = 0; hist < core_data.N_histograms; ++hist) {
        std::cout << ">> Histogram " << hist << ", Property: " << core_data.properties[core_data.output_channels[hist].input_property_index].name
                  << ", Statistic: " << static_cast<int>(core_data.output_channels[hist].statistic)
                  << ", Order: " << core_data.output_channels[hist].order
                  << ", Output Histogram Index: " << core_data.output_channels[hist].output_histogram_index
                  << "\n";
        const int storage_rows = core_data.N_histograms + 2;
        const int n_bins = core_data.output_data.size() / storage_rows;
        for (int bin = 0; bin < std::min(10, n_bins); ++bin) {
            std::cout << core_data.output_data[bin + (hist + 2) * n_bins] << " ";
        }
        std::cout << "\n";
    }

    // Return 0 to indicate successful execution
    return 0;
}