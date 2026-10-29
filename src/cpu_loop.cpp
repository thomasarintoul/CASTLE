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
#include "io.cpp"
#include "data.cpp"
#include "config.cpp"
#include "analysis.cpp"


// defines the loop for CPU computation of the structure function
void process_pairs(std::size_t pair_start, std::size_t pair_end,
                   const std::vector<double>& input_data,
                   std::size_t N_part,
                   const int N_components,
                   const int N_histograms,
                   const std::vector<InputProperty>& properties,
                   std::vector<double>& output_data,
                   std::vector<OutputChannel>& output_channels) {

    // Loop over pairs, calculate i and j from the pair index, and compute the structure function for each pair
    int last_reported_percent = 0;
    for (std::size_t pair_index = pair_start; pair_index < pair_end; ++pair_index) {
        if (pair_index % 10000000 == 0){
            const std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::cout << "Time: " << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S")
                      << " - Processing pair index: " << pair_index << " of " << pair_end << "\n";
        }


        // Calculate i and j from the pair index
        std::size_t i = static_cast<std::size_t>(
            std::floor(
                (2.0 * N_part - 1.0 - std::sqrt((2.0 * N_part - 1.0) * (2.0 * N_part - 1.0) - 8.0 * pair_index)) / 2.0
            )
        );
        std::size_t j = i + 1 + (pair_index  - (i * (2 * N_part - i - 1) / 2));

    //// Loop over all pairs of particles to compute structure functions
    //// for (int i = 0; i < N_part; ++i) {
    ////     for (int j = i + 1; j < N_part; ++j) {

        // Get the positions of particles i and j
        // Note: by definition, x y z are the first three components of the input_data array
        std::array<double, 3> pos_i = {input_data[N_part*0 + i], input_data[N_part*1 + i], input_data[N_part*2 + i]};
        std::array<double, 3> pos_j = {input_data[N_part*0 + j], input_data[N_part*1 + j], input_data[N_part*2 + j]};


        // Compute distance and other relevant quantities between particle i and j
        const int storage_rows = N_histograms + 2;
        const int N_bins = output_data.size() / storage_rows;
        const std::vector<double> bin_edges(output_data.begin(), output_data.begin() + N_bins);

        double distance = calc_separation(pos_i, pos_j);
        int bin_index = find_bin(distance, bin_edges, N_bins); //TODO: define bin_edges and N_bins within this scope
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







void run_cpu_sf(
    const std::vector<double>& input_data,              // 2D array of input data (N_part x N_components)
    std::size_t N_part,                           // number of particles
    const int N_components,                     // number of components in the input data
    const int N_histograms,                     // number of histograms to compute
    const std::vector<InputProperty>& properties,                  // vector of InputProperty structs defining the properties to compute
    std::vector<double>& output_data,                  // 2D array to hold the computed histograms (N_histograms x N_bins)
    std::vector<OutputChannel>& output_channels,  // mapping of output channels to input and output histogram indices
    int N_threads = 1
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

    // TODO: Remove this testing limit later.
    N_part = std::min<std::size_t>(N_part, 100000);


    // Total number of unique pairs (i, j) with i < j: N * (N - 1) / 2.
    // Use 64-bit integers because this overflows int quickly (N = 100,000 gives ~5 billion).
    const std::uint64_t num_particles = static_cast<std::uint64_t>(N_part);
    const std::uint64_t total_pairs   = num_particles * (num_particles - 1) / 2;
    int last_reported_percent = -1;  // -1 so the first iteration always prints

    // If N_threads > 1, we will divide the total pairs into chunks for each thread to process
    std::vector<std::vector<double>> thread_output_data;
    if (N_threads > 1) {
        // Create local output arrays for each thread to avoid race conditions
        // Each local array is initialised to zero
        thread_output_data.assign(N_threads, std::vector<double>(output_data.size(), 0.0));

        // Copy bin edges from the main output_data array to each thread's local output array
        for (int thread_id = 0; thread_id < N_threads; ++thread_id) {
            std::copy(output_data.begin(), output_data.begin() + N_bins, thread_output_data[thread_id].begin());
        }
    } else {
        // If only one thread, we can use the main output_data array directly
        thread_output_data.assign(1, output_data);
    }

    #pragma omp parallel num_threads(N_threads)
    {
        // Get the thread ID and total number of threads
        const int thread_id = omp_get_thread_num();
        const int total_threads = omp_get_num_threads();
        // Calculate number of pairs per thread (excl. remainder)
        const std::size_t pairs_per_thread = total_pairs / total_threads;
        // Calculate the remainder pairs that need to be distributed among the first few threads
        const std::size_t remainder_pairs = total_pairs % total_threads;
        
        const std::size_t pair_start = thread_id * pairs_per_thread + \
                                       std::min<std::size_t>(thread_id, remainder_pairs);
        const std::size_t pair_end = pair_start + pairs_per_thread + \
                                     (thread_id < remainder_pairs ? 1 : 0);

        std::cout << "Thread " << thread_id << " processing pairs from " << pair_start << " to " << pair_end << "\n";

        process_pairs(pair_start, pair_end,
                 input_data,
                 N_part,
                 N_components,
                 N_histograms,
                 properties,
                 thread_output_data[thread_id],  // Use the current thread's output data
                 output_channels);


    }



    // After processing, combine results of all threads into the main output_data array
    // In the case of single threaded, this effectively copies the results from thread_output_data[0] to output_data
    // For multi-threaded, this sums the results from all threads into the main output_data array
    for (int thread_id = 0; thread_id < N_threads; ++thread_id) {
        // Loop over all rows except the bin edges (row 0), which are already handled
        for (std::size_t idx = N_bins; idx < output_data.size(); ++idx) {
            output_data[idx] += thread_output_data[thread_id][idx];
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

    int N_threads = core_data.config_analysis["analysis"]["n_threads"].value_or(1);
    // Call the function to run the CPU structure function computation
    run_cpu_sf(core_data.input_data,         // 2D array of input data (N_part x N_components)
               core_data.N_part,             // N of elements in the input data (number of particles)
               core_data.N_components,       // N of input data components (3 for vectors, 1 for scalars)
               core_data.N_histograms,       // N of output hists (i.e. number stat-order-prop combinations)
               core_data.properties,         // vector of InputProperty structs defining input properties
               core_data.output_data,        // 2D array to hold the computed histograms (N_histograms x N_bins)
               core_data.output_channels,
               N_threads);   // mapping of output channels to input and output histogram indices


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
    for (std::size_t hist = 0; hist < core_data.output_channels.size(); ++hist) {
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

    write_to_files(core_data.output_data, 
                   core_data.output_data.size() / (core_data.N_histograms + 2), 
                   core_data.N_histograms, 
                   core_data.properties, 
                   core_data.output_channels, 
                   core_data.config_output);


    // Return 0 to indicate successful execution
    return 0;
}