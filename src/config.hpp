// TOML Reading
#include <iostream>
#include "toml.hpp"


/**
* Load double-precision values from a binary file.
*
* @param filename:         Path to the binary file.
* @param expected_values:  Expected number of values, or zero to skip validation.
* @return:                 The values read from the file.
* @throws:                 std::runtime_error If the file is invalid, cannot be read, or has an
*                          unexpected number of values.
*/
std::vector<double> load_binary_doubles(const std::string& filename, std::size_t expected_values);

/**
* Read a TOML file and return its contents as a toml::table.
*
* @param filename:         Path to the TOML file.
* @return:                 toml::table containing the parsed TOML data.
* @throws:                 std::runtime_error If the file cannot be opened or parsed.
*/
toml::table read_toml_file(const std::string& filename);

/**
* Calculate the total number of components based on the provided fields.
*
* @param fields:           toml::array containing the field definitions.
* @return:                 Total number of components.
*/
int calc_N_components(toml::array& fields);

/**
* Calculate the total number of histograms based on the provided fields.
*
* @param fields:           toml::array containing the field definitions.
* @return:                 Total number of histograms.
*/
int calc_N_histograms(toml::array& fields);
