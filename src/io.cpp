#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "io.hpp"
#include "toml.hpp"
#include "data.hpp"
// existing load_bin_*, save_to_binary, linspace, logspace functions

namespace {

const char* statistic_name(Statistic statistic) {
    switch (statistic) {
        case Statistic::ScalarDifference: return "ScalarDifference";
        case Statistic::VectorDifference: return "VectorDifference";
        case Statistic::Parallel: return "Parallel";
        case Statistic::Perpendicular: return "Perpendicular";
    }
    return "Unknown";
}

}

void in_out_human_summary(const std::vector<InputProperty>& properties,
                          const std::vector<OutputChannel>& output_channels,
                          const std::vector<double>& output_histograms,
                          int N_components,
                          int N_histograms,
                          int N_bins){
    // Print human-readable summary of the input and output data

    std::cout << "\n========== Data summary ==========" << "\n"
              << "Input data: "<< N_components << " components (row-major)\n";

    // Loop through each property and print its details along with the corresponding output mappings
    for (std::size_t property_index = 0; property_index < properties.size(); ++property_index) {
        const InputProperty& property = properties[property_index];

        std::cout << "Property " << property_index << ": " << property.name
                  << " (" << (property.type == PropertyType::Scalar ? "scalar" : "vector")
                  << ")" << " (file: " << property.filename << ")\n  input_data component columns: ";

        for (std::size_t component : property.components) {
            std::cout << component;
            std::cout << " ";
        }
        std::cout << "\n  requested output mappings:\n";

        bool mapped = false;
        for (const auto& channel : output_channels) {
            if (channel.input_property_index != property_index) {
                continue;
            }
            mapped = true;
            std::cout << "    statistic = " << statistic_name(channel.statistic)
                      << ", order = " << channel.order
                      << " -> output_array[" << N_histograms
                      << " * " << channel.output_histogram_index << " + bin]"
                      << " (histogram channel " << channel.output_histogram_index << ")\n";
        }
        if (!mapped) {
            std::cout << "    none\n";
        }
    }

    std::cout << "Output data: " << N_bins << " bins x " << N_histograms
              << " histogram channels, " << output_histograms.size()
              << " values total\n"
              << "==================================\n";
}

void save_to_binary(const std::string& filename, const std::vector<double>& data) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }
    
    file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(double));
    file.close();
}



void write_to_files(const std::vector<double>& output_data, 
                    int N_bins, int N_histograms, 
                    const std::vector<InputProperty>& properties,
                    const std::vector<OutputChannel>& output_channels,
                    toml::table& config_output) {
    // Write the output data to binary files for each histogram channel
    const toml::table* output_config = config_output["output"].as_table();
    if (!output_config) {
        throw std::runtime_error("Missing [output] configuration table.");
    }

    const auto output_base_dir = output_config->get_as<std::string>("output_base_dir");
    const auto output_sub_dir = output_config->get_as<std::string>("output_sub_dir");
    const auto output_file_prefix = output_config->get_as<std::string>("output_file_prefix");
    if (!output_base_dir || !output_sub_dir || !output_file_prefix) {
        throw std::runtime_error(
            "[output] must define output_base_dir, output_sub_dir, and output_file_prefix.");
    }

    if (N_bins < 2 || N_histograms < 0 ||
        output_data.size() != static_cast<std::size_t>(N_bins) *
                                  static_cast<std::size_t>(N_histograms + 2)) {
        throw std::runtime_error("Output data size does not match the histogram configuration.");
    }

    if (output_channels.size() > static_cast<std::size_t>(N_histograms)) {
        throw std::runtime_error("More output channels were provided than histograms.");
    }

    std::string output_dir = output_base_dir->get();
    output_dir += output_sub_dir->get();

    const std::string filename_prefix = output_file_prefix->get();

    // Start with bin edges and counts - these are always the first two rows in the output_data array
    std::vector<double> bin_edges(output_data.begin(), output_data.begin() + N_bins);
    std::vector<double> counts(output_data.begin() + N_bins, output_data.begin() + 2 * N_bins);

    // Write bin edges and counts to binary files
    save_to_binary(output_dir + filename_prefix + "_bin_edges.bin", bin_edges);
    save_to_binary(output_dir + filename_prefix + "_counts.bin", counts);

    // Loop over the property-statistic-order combinations, which are stored in the subsequent rows of the output_data array
    for (std::size_t hist = 0; hist < output_channels.size(); ++hist) {
        const OutputChannel& channel = output_channels[hist];
        if (channel.input_property_index >= properties.size() ||
            channel.output_histogram_index >= static_cast<std::size_t>(N_histograms)) {
            throw std::runtime_error("Output channel contains an invalid index.");
        }
        const InputProperty& property = properties[channel.input_property_index];

        const std::size_t histogram_index = channel.output_histogram_index;
        std::vector<double> histogram(
            output_data.begin() + (histogram_index + 2) * N_bins,
            output_data.begin() + (histogram_index + 3) * N_bins);
        save_to_binary(output_dir + \
                       filename_prefix + "_" + property.name + "_" +\
                       statistic_name(channel.statistic) +\
                       "_order_" + std::to_string(channel.order) + ".bin",\
                       histogram);
        }
}