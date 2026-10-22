#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "io.hpp"
#include "toml.hpp"
#include "data.hpp"
// existing load_bin_*, save_to_binary, linspace, logspace functions



void in_out_human_summary(const std::vector<InputProperty>& properties,
                          const std::vector<OutputChannel>& output_channels,
                          const std::vector<double>& output_histograms,
                          int N_components,
                          int N_histograms,
                          int N_bins){
    // Print human-readable summary of the input and output data

    /// switch to convert Statistic enum to string for printing
    const auto statistic_name = [](Statistic statistic) {
        switch (statistic) {
            case Statistic::ScalarDifference: return "ScalarDifference";
            case Statistic::VectorDifference: return "VectorDifference";
            case Statistic::Parallel: return "Parallel";
            case Statistic::Perpendicular: return "Perpendicular";
        }
        return "Unknown";
    };


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