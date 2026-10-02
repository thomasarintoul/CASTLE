// TOML Reading
#include <iostream>
#include "toml.hpp"

toml::table read_toml_file(const std::string& filename);

int calc_N_components(toml::array& fields);
int calc_N_histograms(toml::array& fields);
