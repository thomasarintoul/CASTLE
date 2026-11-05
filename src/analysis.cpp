#include "compatibility.h"
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <vector>

HD double calc_separation(const std::array<double, 3>& pos_i, const std::array<double, 3>& pos_j) {
    /**
     * Calculate the Euclidean distance between two 3D positions.
     *
     * @param pos_i: 3D position of particle i.
     * @param pos_j: 3D position of particle j.
     * @return: Euclidean distance between pos_i and pos_j.
     */
    double dx = pos_i[0] - pos_j[0];
    double dy = pos_i[1] - pos_j[1];
    double dz = pos_i[2] - pos_j[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

HD std::array<double, 3> get_position_vector(const double* input_data, 
                                        const size_t N_part, 
                                        const size_t i) {
    /**
     * Extract the position vector for a given particle index.
     *
     * @param input_data    : 1D vector containing the input data (flattened 2D array).
     * @param N_part        : Total number of particles.
     * @param i             : Index of the particle.
     * @return              : A vector containing the position components [x, y, z] for particle i.
     */

    // We assume that the first three components of input data are
    // position by construction in config.cpp
    return {input_data[N_part * 0 + i], input_data[N_part * 1 + i], input_data[N_part * 2 + i]};
}

HD std::array<double, 3> calc_sep_vector(const double* input_data, 
                                    const size_t N_part, 
                                    const size_t i, const size_t j) {
    /**
     * Calculate the separation vector between two particles.
     *
     * @param input_data    : 1D vector containing the input data (flattened 2D array).
     * @param N_part        : Total number of particles.
     * @param i             : Index of the first element.
     * @param j             : Index of the second element.
     * @return              : A vector containing the separation vector components [dx, dy, dz].
     */

     // We assume that the first three components of input data are
     // position by construction in config.cpp
    double dx = input_data[N_part * 0 + i] - input_data[N_part * 0 + j];
    double dy = input_data[N_part * 1 + i] - input_data[N_part * 1 + j];
    double dz = input_data[N_part * 2 + i] - input_data[N_part * 2 + j];

    return {dx, dy, dz};
}

HD int find_bin(double separation, const double* bins, int n_bins) {

    if (separation < bins[0] || separation >= bins[n_bins - 1]) return -1;

    int left = 0;
    int right = n_bins - 1;

    while (right - left > 1) {
        int mid = (left + right) / 2;

        if (separation >= bins[mid])
            left = mid;
        else
            right = mid;
    }

    return left;
}

HD double compute_statistic(const double* input_data,
                         const size_t N_part,
                         const size_t i,
                         const size_t j,
                         const std::size_t* input_components,
                         const Statistic statistic,
                         const int order) {
    /**
     * Compute the specified statistic for a given pair of particles and a specific property.
     * @param input_data            : 1D vector containing the input data (flattened 2D array).
     * @param N_part                : Total number of particles.
     * @param i                     : Index of the first particle.
     * @param j                     : Index of the second particle.
     * @param input_components      : Columns of the input data.
     * @param statistic             : The statistic to compute (e.g., ScalarDifference).
     * @param order                 : The order of the statistic (e.g., 2 for second order, 3 for third order).
     * @return                      : The computed statistic value for the given pair of particles and property.
     */

    // Select appropriate computation based on the statistic type
    switch (statistic) {
        case Statistic::ScalarDifference: {
            // For scalar properties, compute the absolute difference raised to the specified order
            double value_i = input_data[N_part * input_components[0] + i];
            double value_j = input_data[N_part * input_components[0] + j];

            double difference = std::abs(value_i - value_j);
            return std::pow(difference, order);
        }
        case Statistic::VectorDifference: {
            // For vector properties, compute the Euclidean distance between the two vectors raised to the specified order
            double value_i_x = input_data[N_part * input_components[0] + i];
            double value_i_y = input_data[N_part * input_components[1] + i];
            double value_i_z = input_data[N_part * input_components[2] + i];

            double value_j_x = input_data[N_part * input_components[0] + j];
            double value_j_y = input_data[N_part * input_components[1] + j];
            double value_j_z = input_data[N_part * input_components[2] + j];

            double dx = value_i_x - value_j_x;
            double dy = value_i_y - value_j_y;
            double dz = value_i_z - value_j_z;

            return std::pow(std::sqrt(dx * dx + dy * dy + dz * dz), order);
        }
        case Statistic::Parallel: {
            // For parallel statistics, compute the component of the difference along the separation vector raised to the specified order
            double value_i_x = input_data[N_part * input_components[0] + i];
            double value_i_y = input_data[N_part * input_components[1] + i];
            double value_i_z = input_data[N_part * input_components[2] + i];

            double value_j_x = input_data[N_part * input_components[0] + j];
            double value_j_y = input_data[N_part * input_components[1] + j];
            double value_j_z = input_data[N_part * input_components[2] + j];

            // Compute the separation vector
            std::array<double, 3> sep_vector = calc_sep_vector(input_data, N_part, i, j);

            // Compute the dot product of the difference vector and the separation vector
            double dot_product = (value_i_x - value_j_x) * sep_vector[0] +
                                 (value_i_y - value_j_y) * sep_vector[1] +
                                 (value_i_z - value_j_z) * sep_vector[2];
            // Compute the magnitude of the separation vector
            double sep_magnitude = std::sqrt(sep_vector[0] * sep_vector[0] + 
                                             sep_vector[1] * sep_vector[1] + 
                                             sep_vector[2] * sep_vector[2]);

            // Compute the parallel component
            double parallel_component = dot_product / sep_magnitude;
            return std::pow(parallel_component, order);
        }
        case Statistic::Perpendicular: {
            // For perpendicular statistics, compute the component of the difference perpendicular to the separation vector raised to the specified order
            double value_i_x = input_data[N_part * input_components[0] + i];
            double value_i_y = input_data[N_part * input_components[1] + i];
            double value_i_z = input_data[N_part * input_components[2] + i];

            double value_j_x = input_data[N_part * input_components[0] + j];
            double value_j_y = input_data[N_part * input_components[1] + j];
            double value_j_z = input_data[N_part * input_components[2] + j];

            // Compute the separation vector
            std::array<double, 3> sep_vector = calc_sep_vector(input_data, N_part, i, j);
        
            // Compute the difference vector
            std::array<double, 3> diff_vector = {value_i_x - value_j_x, value_i_y - value_j_y, value_i_z - value_j_z};

            // Get a vector perpendicular to the separation vector (in 2D, we can use a simple rotation)
            std::array<double, 3> normal_vector = {-sep_vector[1], sep_vector[0], 0.0};

            // Compute perpendicular component
            double perpendicular_component = std::sqrt(
                (diff_vector[0] * normal_vector[0]) * (diff_vector[0] * normal_vector[0]) +
                (diff_vector[1] * normal_vector[1]) * (diff_vector[1] * normal_vector[1]) +
                (diff_vector[2] * normal_vector[2]) * (diff_vector[2] * normal_vector[2])
            );
            return std::pow(perpendicular_component, order);
        }
        default:
            return 0.0;
    }
}


HD Pair_Index calc_i_j(std::size_t pair_index, std::size_t N_part) {
    /**
     * Calculate the particle indices i and j from a given pair index.
     *
     * @param pair_index: The index of the particle pair.
     * @param N_part: Total number of particles.
     * @return: A vector containing the indices [i, j] of the particle pair.
     */
    std::size_t i = static_cast<std::size_t>(
        std::floor(
            (2.0 * N_part - 1.0 - std::sqrt((2.0 * N_part - 1.0) * (2.0 * N_part - 1.0) - 8.0 * pair_index)) / 2.0
        )
    );
    std::size_t j = i + 1 + (pair_index  - (i * (2 * N_part - i - 1) / 2));
    return {i, j};
}

HD std::array<double, 3> get_position(const double* input_data, 
                                    const size_t N_part, 
                                    const size_t i) {
    /**
     * Extract the position vector for a given particle index.
     *
     * @param input_data    : 1D vector containing the input data (flattened 2D array).
     * @param N_part        : Total number of particles.
     * @param i             : Index of the particle.
     * @return              : An array containing the position components [x, y, z] for particle i.
     */

    // We assume that the first three components of input data are
    // position by construction in config.cpp
    return {input_data[N_part * 0 + i], input_data[N_part * 1 + i], input_data[N_part * 2 + i]};
}

H void process_pairs(std::size_t pair_start, std::size_t pair_end,
                   const double* input_data,
                   std::size_t N_part,
                   const int N_components,
                   const int N_histograms,
                   const std::vector<InputProperty>& properties,
                   std::vector<double>& output_data,
                   std::vector<OutputChannel>& output_channels) {

    // Loop over pairs, calculate i and j from the pair index, and compute the structure function for each pair
    // int last_reported_percent = 0;
    
    // Compute distance and other relevant quantities between particle i and j
    const int storage_rows = N_histograms + 2;
    const int N_bins = output_data.size() / storage_rows;
    const double* bin_edges = output_data.data();

    for (std::size_t pair_index = pair_start; pair_index < pair_end; ++pair_index) {
        if (pair_index % 10000000 == 0){
            const std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::cout << "Time: " << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S")
                      << " - Processing pair index: " << pair_index << " of " << pair_end << "\n";
        }


        // Calculate i and j from the pair index
        Pair_Index indices = calc_i_j(pair_index, N_part);
        std::size_t i = indices.i;
        std::size_t j = indices.j;

        if (i >= N_part || j >= N_part || j <= i) {
            continue;
        }

        // Get the positions of particles i and j
        // Note: by definition, x y z are the first three components of the input_data array
        std::array<double, 3> pos_i = get_position(input_data, N_part, i);
        std::array<double, 3> pos_j = get_position(input_data, N_part, j);



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
            // Get config values for prop-stat-order calculation
            std::size_t input_property_index = channel.input_property_index;
            Statistic statistic = channel.statistic;
            int order = channel.order;
            std::size_t output_histogram_index = channel.output_histogram_index;
            const std::size_t* input_components = properties[input_property_index].components.data();

            // Compute the relevant quantity based on the statistic and order
            double value = compute_statistic(input_data, N_part, i, j, input_components, statistic, order);
            
            // Add the computed value to the appropriate histogram in the output_data array
            output_data[N_bins * (2 + output_histogram_index) + bin_index] += value;
        }
    }
}
