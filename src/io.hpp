#pragma once

#include <vector>

#include "data.hpp"

// existing load_bin_*, save_to_binary, linspace, logspace functions

void in_out_human_summary(const std::vector<InputProperty>& properties,
                          const std::vector<OutputChannel>& output_channels,
                          const std::vector<double>& output_histograms,
                          int N_components,
                          int N_histograms,
                          int N_bins);

void save_to_binary(const std::string& filename, const std::vector<double>& data);

void write_to_files(const std::vector<double>& output_data, 
                    int N_bins, int N_histograms, 
                    const std::vector<InputProperty>& properties,
                    const std::vector<OutputChannel>& output_channels,
                    toml::table& config_output);