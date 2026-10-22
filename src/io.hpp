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