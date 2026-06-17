#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <iterator>
#include <numeric>
#include <tuple>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <atomic>
#include <array>
#include <thread>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <array>
#include <unordered_map>
// #include "cnpy.h"
// #include <Eigen/Dense> // For vector operations
#include </home/user/c23052501/local/include/tqdm.h> // progress bar

using namespace std;


std::vector<double> load_bin(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary | std::ios::ate); // Open at end
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    size_t num_bytes = file.tellg();  // Get total file size
    file.seekg(0, std::ios::beg);  // Move back to the beginning

    if (num_bytes % (sizeof(double) * 3) != 0) {
        throw std::runtime_error("File size is not a multiple of 3 * sizeof(double)");
    }

    size_t num_elements = num_bytes / sizeof(double);  // Total doubles
    std::vector<double> data(num_elements);

    file.read(reinterpret_cast<char*>(data.data()), num_bytes);
    return data;
}


std::vector<std::array<double, 3>> load_bin_3d(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    size_t num_bytes = file.tellg();
    file.seekg(0, std::ios::beg);

    if (num_bytes % (sizeof(double) * 3) != 0) {
        throw std::runtime_error("File size is not a multiple of 3 * sizeof(double)");
    }

    size_t N = num_bytes / (sizeof(double) * 3);
    std::vector<std::array<double, 3>> data(N);

    file.read(reinterpret_cast<char*>(data.data()), num_bytes);
    return data;
}


vector<double> linspace(double start, double end, int num) {
    vector<double> linspaced;
    double delta = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        linspaced.push_back(start + delta * i);
    }
    return linspaced;
}


vector<double> logspace(double start, double end, int num) {
    vector<double> logspaced;
    double delta = (log10(end) - log10(start)) / (num - 1);
    for (int i = 0; i < num; ++i) {
        logspaced.push_back(pow(10, log10(start) + delta * i));
    }
    return logspaced;
}

double three_dim_sep(std::array<double, 3> vec1, std::array<double, 3> vec2){
    double sum = 0.0;
    for (int i = 0; i < 3; i++){
        sum += pow(vec1[i] - vec2[i], 2);
    }
    return sqrt(sum);
}

pair<vector<double>, vector<double>> VSF(std::vector<std::array<double, 3>> gas_positions, 
                                         std::vector<std::array<double, 3>> gas_velocities, 
                                         int n_bins, double bin_min, double bin_max, 
                                         bool log_space = true, bool norm = true, 
                                         bool subset = true, int subset_size=1000,
                                         int n_threads=1){
    std::cout << "Computing VSF..." << std::endl;

    unsigned long long N = gas_positions.size();
    std::cout << "Number of particles: " << N << std::endl;

    if(subset){
        std::cout << "Subsetting to " << subset_size << " particles." << std::endl;
        // equiv of np.arange()
        vector<int> indices(N);
        iota(indices.begin(), indices.end(), 0);
        // random number generator
        random_device rd;
        mt19937 g(rd());
        // shuffle indices into a random order
        shuffle(indices.begin(), indices.end(), g);
        // select the first subset_size indices
        indices.resize(subset_size);

        std::vector<std::array<double, 3>> selected_positions(subset_size);
        std::vector<std::array<double, 3>> selected_velocities(subset_size);
        for (int i = 0; i < subset_size; ++i) {
        selected_positions[i] = gas_positions[indices[i]];
        selected_velocities[i] = gas_velocities[indices[i]];
        }
        gas_positions = selected_positions;
        gas_velocities = selected_velocities;
        N = subset_size;
    }

    vector<double> bins;
    if (log_space){
        bins = logspace(bin_min, bin_max, n_bins);
    }
    else{
        bins = linspace(bin_min, bin_max, n_bins);
    }

    vector<double> bin_hist(n_bins, 0.0);
    vector<double> VSF_data(n_bins, 0.0);

    if (n_threads > 1){
        omp_set_num_threads(n_threads);
    }

    // Initialise tqdm progress bar
    tqdm bar;
    bar.set_theme_line();
    bar.set_label("Computing VSF");
    unsigned long long total_iterations = N * (N - 1) / 2;
    std::cout << "Total iterations: " << total_iterations << std::endl;
    unsigned long long current_iteration = 0;

    //TODO: Flatten this loop as in the GPU version
    if (n_threads > 1){
    #pragma omp parallel for schedule(dynamic) reduction(+:current_iteration)
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double r = three_dim_sep(gas_positions[i], gas_positions[j]);
            double v = pow(three_dim_sep(gas_velocities[i], gas_velocities[j]), 2);
            auto it = lower_bound(bins.begin(), bins.end(), r);
            int index = distance(bins.begin(), it);
            if (index < n_bins) {
                #pragma omp atomic
                bin_hist[index] += 1;
                #pragma omp atomic
                VSF_data[index] += v;
            }
            #pragma omp atomic
            current_iteration++;

        }
    }
    }
    else{
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                double r = three_dim_sep(gas_positions[i], gas_positions[j]);
                double v = pow(three_dim_sep(gas_velocities[i], gas_velocities[j]), 2);
                auto it = lower_bound(bins.begin(), bins.end(), r);
                int index = distance(bins.begin(), it);
                if (index < n_bins) {
                    bin_hist[index] += 1;
                    VSF_data[index] += v;
                }
                current_iteration++;
                bar.progress_huge(current_iteration, total_iterations);
            }
        }
        bar.finish();
    }

    if (norm){
        for (int i = 0; i < n_bins; i++){
            if (bin_hist[i] > 0){
                VSF_data[i] /= bin_hist[i];
            }
        }
    }

    return make_pair(bins, VSF_data);
    }



void save_to_binary(const std::string& filename, const std::vector<double>& data) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }
    
    file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(double));
    file.close();
}

unordered_map<string, string> read_parameters(const string& param_file) {
    // Create the unordered_map to store parameters (similar to a dict in python)
    unordered_map<string, string> params;
    
    // Open the parameter file
    ifstream file(param_file);
    if (!file.is_open()) {
        cerr << "Error: Unable to open parameter file: " << param_file << endl;
        exit(1);
    }

    // Read the file line by line
    string line;
    while (getline(file, line)) {
        // Skip lines starting with '%'
        if (!line.empty() && line[0] == '%') {
            continue;
        }

        istringstream iss(line);
        string key, value;
        if (iss >> key >> value) {
            params[key] = value;
        }
    }
    file.close();
    return params;
}


int get_subset_size(unordered_map<string, string> params, int N_cells){
    // Check if the subset size is specified in the parameters
    if (params.find("subset_size") != params.end()) {
        // If specified, convert it to an integer
        int subset_size = stoi(params["subset_size"]);
        return subset_size;
    }
    // If subset_size is not specified, check for defined proportions
    else if ((params.find("subset_proportion") != params.end())
            && (params.find("total_proportion") == params.end())) {
            // If subset_proportion is specified, get subset size from this
            double subset_proportion = stod(params["subset_proportion"]);
            int subset_size = static_cast<int>(N_cells * subset_proportion);
            
            return subset_size;
            }
    else if ((params.find("total_proportion") != params.end())
             && (params.find("subset_proportion") == params.end())) {
                // If total_proportion is specified, provided
                // Check we have provided a number of samples
                // to take, then calculate subset size from that
                if (params.find("n_samples") != params.end()){
                    int subset_size = static_cast<int>(N_cells * 
                    (stod(params["total_proportion"]) / 
                    stoi(params["n_samples"])));
                    return subset_size;
                }           
                else {
                    throw std::runtime_error("Error: 'n_samples' is missing, "
                                        "but 'total_proportion' is provided. "
                                        "Both are required together.");
                }
             }    
    else if ((params.find("total_proportion") != params.end())
                && (params.find("subset_proportion") != params.end())) {
        throw std::runtime_error("Error: Both 'total_proportion' and "
                    "'subset_proportion' are provided. Only one should "
                    "be specified if not providing 'subset_size' directly.");
    }
    
    else if ((params.find("total_proportion") == params.end())
                && (params.find("subset_proportion") == params.end())) {
        throw std::runtime_error("Error: Neither 'subset_size' nor "
                    "'total_proportion' nor 'subset_proportion' is "
                    "provided. One of them must be specified.");
    }
    
    else {
        throw std::runtime_error("Error: Unhandled case encountered. "
                                "Please check the input parameters.");
    }
}

