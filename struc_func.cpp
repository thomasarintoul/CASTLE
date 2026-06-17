#include <iostream>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <fstream>
#include "struc_func_functions.hpp"



using namespace std;


int main(int argc, char* argv[]){

    // Raise an error if not enough arguments are provided
    if (argc != 4){
        std::cerr << "Usage: " << argv[0] << " <dataset1 path - e.g. pos> <dataset2 path - e.g. vel> <parameter file>\n";
        return 1;
    }

    std::string pos_path = argv[1];
    std::string vel_path = argv[2];

    // Load in position and velocity data
    std::vector<std::array<double, 3>> positions = load_bin_3d(pos_path);
    std::vector<std::array<double, 3>> velocities = load_bin_3d(vel_path);
    // Confirm number of cells loaded
    std::cout << "Detected " << positions.size() << " position rows.\n";
    std::cout << "Detected " << velocities.size() << " velocity rows.\n";

    // Read parameters from the parameter file
    std::string param_file = argv[3];
    unordered_map<string, string> params = read_parameters(param_file);
    

    int n_samples = 1;
    if (params.find("n_samples") != params.end()){
        n_samples = stoi(params["n_samples"]);
    }
    
    int n_bins = stoi(params["n_bins"]);
    double bin_min = stod(params["bin_min"]);
    double bin_max = stod(params["bin_max"]);
    bool log_space = (params["log_space"] == "true");
    bool norm = (params["norm"] == "true");
    bool subset = (params["subset"] == "true");
    int n_threads = stoi(params["n_threads"]);
    std::cout << "loaded params.\n";
    
    
    int subset_size = get_subset_size(params, positions.size());

    // Print out loaded parameters
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "-------   LOADED PARAMETERS   -------" << std::endl;
    std::cout << " Note: for boolean values, 1 == True " << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "- n_samples:                 " << n_samples << std::endl;
    std::cout << "- n_bins:                    " << n_bins << std::endl;
    std::cout << "- bin_min:                   " << bin_min << std::endl;
    std::cout << "- bin_max:                   " << bin_max << std::endl;
    std::cout << "- log_space:                 " << log_space << std::endl;
    std::cout << "- norm:                      " << norm << std::endl;
    std::cout << "- subset:                    " << subset << std::endl;
    std::cout << "- subset_size:               " << subset_size << std::endl;
    std::cout << "- subset_size * n_samples:   " << subset_size * n_samples << std::endl;
    std::cout << "- n_threads:                 " << n_threads << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    
    for (int i = 0; i < n_samples; i++){
        std::cout << "Sample " << i << ". Total samples:" << n_samples << std::endl;
        // Compute VSF
        auto result = VSF(positions, velocities, 
            n_bins, bin_min, bin_max,  // n_bins, bin_min, bin_max
            log_space, norm,    // log_space, norm
            subset, subset_size,   // subset, subset_size
            n_threads);            // n_threads

        // Separate bins and VSF data
        vector<double> bins = result.first;
        vector<double> VSF_data = result.second;

        std::string iter = std::to_string(i);
        save_to_binary("tmp_files/vsf_out/"+
            params["output_file_prefix"]+
            "_bin_data_s"+iter+".bin", bins);
        save_to_binary("tmp_files/vsf_out/"+
            params["output_file_prefix"]+
            "_VSF_data_s"+iter+".bin", VSF_data);
    }

    return 0;
}

