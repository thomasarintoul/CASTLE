// The original main file but shrunk considerably
// It now just reads the input files and calls the CPU or CUDA computation function.

#include <iostream>
#include <vector>
#include <array>
#include <unordered_map>
#include <cuda_runtime.h>
#include "struc_func_functions.cuh" 


#if defined(__CUDACC__)
    #define SF_HD __host__ __device__
#else
    #define SF_HD
#endif

int main(int argc, char* argv[]) {
        
    // Call functions from config.cpp to read the TOML files for 
    // analysis parameters, output parameters,
    // and to read the property fields in use
    // Gives us:
    // - a set of analysis parameters (e.g., number of bins, 
    //   min/max separation, etc.)
    // - a set of output parameters (e.g., output file names, etc.)
    // - a set of property fields and observations to loop over



    // Call io.cpp files to read in the input data to large array



    // Compile all relevant information into a "Problem" 
    // struct (rename this) - defined in pair_kernel.hpp



    // Pass Problem struct to applicable compute function (CPU or 
    // CUDA) - these are in cpu_loop.cpp and cuda_loop.cu



    // Unpack returned histograms into relevant output files (e.g., pair 
    // counts, structure functions, etc.) - defined in pair_kernel.hpp 
    // and io.cpp



    // Clean up memory and exit


    

    return 0;
}