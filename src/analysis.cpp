double calc_separation(const std::array<double, 3>& pos_i, const std::array<double, 3>& pos_j) {
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

std::vector<double> get_position_vector(const std::vector<double>& input_data, 
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


std::vector<double> calc_sep_vector(const std::vector<double>& input_data, 
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



int find_bin(double separation, const std::vector<double>& bins, int n_bins) {

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



double compute_statistic(const std::vector<double>& input_data,
                         const size_t N_part,
                         const size_t i,
                         const size_t j,
                         const std::vector<size_t> input_components,
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
            std::vector<double> sep_vector = calc_sep_vector(input_data, N_part, i, j);

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
            std::vector<double> sep_vector = calc_sep_vector(input_data, N_part, i, j);
        
            // Compute the difference vector
            std::vector<double> diff_vector = {value_i_x - value_j_x, value_i_y - value_j_y, value_i_z - value_j_z};

            // Get a vector perpendicular to the separation vector (in 2D, we can use a simple rotation)
            std::vector<double> normal_vector = {-sep_vector[1], sep_vector[0], 0.0};

            // Compute perpendicular component
            double perpendicular_component = std::sqrt(
                (diff_vector[0] * normal_vector[0]) * (diff_vector[0] * normal_vector[0]) +
                (diff_vector[1] * normal_vector[1]) * (diff_vector[1] * normal_vector[1]) +
                (diff_vector[2] * normal_vector[2]) * (diff_vector[2] * normal_vector[2])
            );
            return std::pow(perpendicular_component, order);
        }
        default:
            throw std::invalid_argument("Unknown statistic type.");
    }
}
