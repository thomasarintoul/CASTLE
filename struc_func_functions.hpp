#ifndef VSF_FUNCTIONS_H
    #define VSF_FUNCTIONS_H

    #include <iostream>     // strings
    #include <vector>       // vectors
    #include <algorithm>    // vector ops on arrays
    #include <unordered_map> // hash maps

    // #include <vector>
    // #include <cmath>
    // #include <random>
    // #include <iterator>
    // #include <numeric>
    // #include <tuple>
    // #include <fstream>
    // #include "cnpy.h"
    // #include <Eigen/Dense> // For vector operations
    // #include </home/user/c23052501/local/include/tqdm.h> // progress bar

    using namespace std;


    std::vector<double> load_bin(const std::string& filename);


    std::vector<std::array<double, 3>> load_bin_3d(const std::string& filename);


    vector<double> linspace(double start, double end, int num);


    vector<double> logspace(double start, double end, int num);

    double three_dim_sep(std::array<double, 3> vec1, std::array<double, 3> vec2);

    pair<vector<double>, vector<double>> VSF(std::vector<std::array<double, 3>> gas_positions, 
                                            std::vector<std::array<double, 3>> gas_velocities, 
                                            int n_bins, double bin_min, double bin_max, 
                                            bool log_space = true, bool norm = true, 
                                            bool subset = true, int subset_size=1000,
                                            int n_threads=1);

    void save_to_binary(const std::string& filename, const std::vector<double>& data);

    unordered_map<string, string> read_parameters(const string& param_file);

    int get_subset_size(unordered_map<string, string> params, int N_cells);

#endif // VSF_FUNCTIONS_H