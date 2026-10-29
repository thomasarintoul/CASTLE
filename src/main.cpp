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
#include "cpu_loop.cpp"


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

    // TODO: Select backend (CPU or GPU [initially CUDA only but later for AMD/Intel])

    
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