#pragma once

#include <algorithm>  // std::find_if
#include <cstddef>    // std::size_t
#include <cstdint>    // std::int64_t
#include <iostream>   // std::cerr
#include <iterator>   // std::distance
#include <stdexcept>  // std::runtime_error
#include <string>     // std::string
#include <vector>     // std::vector

#include "toml.hpp"



enum class Statistic {
    ScalarDifference,   // Used for 1D scalar properties (e.g., temperature, density, |B|)
    // Note: When the python frontend is implemented, we will need to enable it to create e.g. |B| if B is provided for ScalarDifference
    VectorDifference,   // Used for vector properties (e.g., velocity, B field)
    Parallel,           // Used for vector properties, computes the component of the difference along the separation vector
    Perpendicular       // Used for vector properties, computes the component of the difference perpendicular to the separation vector
};

enum class PropertyType {
    Scalar, // 1 component
    Vector  // 3 components
};

struct InputProperty {
    std::string name;               // Name of the property (e.g., "velocity", "temperature")
    PropertyType type;              // Type of the property (e.g., "scalar", "vector")
    std::string filename;           // Binary file containing the property values
    std::vector<std::size_t> components;                 // Component "columns" in the input array (1 for scalar, 3 for vector)
};

struct OutputChannel{
    std::size_t input_property_index;  // Index of the input property in the input array
    Statistic statistic;                // Statistic to compute for this property (e.g., ScalarDifference, VectorDifference, Parallel, Perpendicular)
    int order;                          // Order of the statistic (e.g., 2 for second order, 3 for third order)
    std::size_t output_histogram_index; // Index of the output histogram "column" in the output array
};

inline void Properties(std::vector<InputProperty>& properties, const toml::array& fields);

inline std::vector<OutputChannel> map_output_channels(const std::vector<InputProperty>& properties, const toml::array& fields);