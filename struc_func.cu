#include <iostream>
#include <vector>
#include <array>
#include <unordered_map>
#include <cuda_runtime.h>
#include "struc_func_functions.cuh"

// Compilation Command on Anubis (Nvidia Titan Xp): nvcc -arch=sm_61 struc_func_mGPU.cu struc_func_functions_mGPU.cu -o vsc_cuda.out

// Full compilation command for multiple files:
// nvcc -O3 \
//   -gencode arch=compute_61,code=sm_61 \
//   -gencode arch=compute_89,code=sm_89 \
//   -gencode arch=compute_90,code=sm_90 \
//   -gencode arch=compute_90,code=compute_90 \   # PTX fallback for anything newer than H100
//   struc_func_mGPU.cu struc_func_functions_mGPU.cu -o vsc_cuda.out

int main(int argc, char* argv[]){

    // Command line arguments
    // -pos / -p "path to position data"
    // -vel / -v "path to velocity data"
    // -scalar / -s "path to scalar data"
    // -an_param -ap "path to analysis parameter file"
    // -io_param -ip "path to I/O parameter file"
    // -snap "snapshot number (integer)"
    // -prop "scalar property to calculate structure function for (e.g. 'velocity', 'temperature')"

    // Default values (empty or sensible defaults)
    std::string pos_path = "";
    std::string an_path  = "";
    std::string io_path  = "";
    int snap             = -1;
    std::string vel_path = "";
    std::string sc_path  = "";
    std::string sc_prop  = "scalar"; // default to velocity structure function if not specified

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-pos" || arg == "-p") && i + 1 < argc) {
            pos_path = argv[++i];
        } else if ((arg == "-vel" || arg == "-v") && i + 1 < argc) {
            vel_path = argv[++i];
        } else if ((arg == "-scalar" || arg == "-s") && i + 1 < argc) {
            sc_path = argv[++i];
        } else if ((arg == "-an_param" || arg == "-ap") && i + 1 < argc) {
            an_path = argv[++i];
        } else if ((arg == "-io_param" || arg == "-ip") && i + 1 < argc) {
            io_path = argv[++i];
        } else if (arg == "-snap" && i + 1 < argc) {
            snap = std::stoi(argv[++i]);
        } else if ((arg == "-prop" || arg == "-P") && i + 1 < argc) {
            sc_prop = argv[++i];
        } else {
            std::cerr << "Unknown or incomplete argument: " << arg << "\n";
            return 1;
        }
    }

    // --------------------------------------------
    // Load in position and velocity data (CPU)
    // --------------------------------------------
    std::vector<std::array<double, 3>> positions = (pos_path.empty()) ? std::vector<std::array<double, 3>>() : load_bin_3d(pos_path);
    std::vector<std::array<double, 3>> velocities = (vel_path.empty()) ? std::vector<std::array<double, 3>>() : load_bin_3d(vel_path);
    std::vector<double> scalars = (sc_path.empty()) ? std::vector<double>() : load_bin_1d(sc_path);

    // --------------------------------------------
    // Confirm number of cells loaded
    // --------------------------------------------
    std::cout << "Detected " << positions.size() << " position rows.\n";
    std::cout << "Detected " << velocities.size() << " velocity rows.\n";
    std::cout << "Detected " << scalars.size() << " scalar rows.\n";
    
    // --------------------------------------------
    // Read parameters from the parameter file
    // --------------------------------------------
    std::string analysis_param_file = (an_path.empty()) ? argv[4] : an_path;
    std::string io_param_file = (io_path.empty()) ? argv[5] : io_path;
    std::unordered_map<std::string, std::string> params = read_parameters(analysis_param_file);
    std::unordered_map<std::string, std::string> io_params = read_parameters(io_param_file);
    
    
  
    int n_bins = stoi(params["n_bins"]);
    double bin_min = stod(params["bin_min"]);
    double bin_max = stod(params["bin_max"]);
    bool log_space = (params["log_space"] == "true");
    bool norm = (params["norm"] == "true");
    
    // Optional toggles
    bool do_vsf = (params["do_vsf"] == "true");
    bool do_vsf_parr = (params["do_vsf_parr"] == "true");
    bool do_vsf_perp = (params["do_vsf_perp"] == "true");
    bool do_scalar = (params["do_scalar"] == "true");
    std::cout << "loaded params.\n";

    // Print out loaded parameters
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "-------   LOADED PARAMETERS   -------" << std::endl;
    std::cout << " Note: for boolean values, 1 == True " << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "- n_bins:                    " << n_bins << std::endl;
    std::cout << "- bin_min:                   " << bin_min << std::endl;
    std::cout << "- bin_max:                   " << bin_max << std::endl;
    std::cout << "- log_space:                 " << log_space << std::endl;
    std::cout << "- norm:                      " << norm << std::endl;
    std::cout << "- do_vsf:                    " << do_vsf << std::endl;
    std::cout << "- do_vsf_parr:               " << do_vsf_parr << std::endl;
    std::cout << "- do_vsf_perp:               " << do_vsf_perp << std::endl;
    std::cout << "- do_scalar:                 " << do_scalar << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    
    // --------------------------------------------
    // Create bins for histogram (CPU)
    // --------------------------------------------
    std::vector<double> bins;
    if (log_space) {
        std::cout << "Creating logarithmic bins...\n";
        // For logspace, we need to ensure that bin_min and bin_max are positive
        if (bin_min <= 0 || bin_max <= 0) {
            std::cerr << "Error: For logarithmic bins, bin_min and bin_max must be positive.\n";
            return 1;
        }
        bins = logspace(bin_min, bin_max, n_bins);
    } else {
        std::cout << "Creating linear bins...\n";
        bins = linspace(bin_min, bin_max, n_bins);
    }

    // --------------------------------------------    
    // Calculate total number of pairs to process
    // --------------------------------------------    
    int N = positions.size();
    long long total_iterations = (long long)N * (N - 1) / 2;
    std::cout << "Total iterations: " << total_iterations << std::endl;



    
    // Convert std::array -> double3
    std::vector<double3> positions_double3(N);
    for (int i = 0; i < N; ++i) {
        positions_double3[i] = make_double3(positions[i][0], positions[i][1], positions[i][2]);
    }
    std::vector<double3> velocities_double3;
    if (!velocities.empty()) {
        velocities_double3.resize(N);
        for (int i = 0; i < N; ++i)
            velocities_double3[i] = make_double3(velocities[i][0], velocities[i][1], velocities[i][2]);
    }

    
    // --------------------------------------------
    // Launch the kernel for N GPUs
    // --------------------------------------------
    int num_GPUs;
    cudaGetDeviceCount(&num_GPUs);
    num_GPUs = std::min(num_GPUs, 8); // Limit to 4 GPUs for this example
    std::cout << "Detected " << num_GPUs << " GPUs. Launching kernel.\n";

    long long pairs_per_gpu = (total_iterations + num_GPUs - 1) / num_GPUs;

    // Per-GPU result buffers on host - these get summed at the end
    std::vector<std::vector<double>> h_bin_hist(num_GPUs, std::vector<double>(n_bins, 0));
    std::vector<std::vector<double>> h_vsf     (num_GPUs, std::vector<double>(n_bins, 0));
    std::vector<std::vector<double>> h_vsf_parr(num_GPUs, std::vector<double>(n_bins, 0));
    std::vector<std::vector<double>> h_vsf_perp(num_GPUs, std::vector<double>(n_bins, 0));
    std::vector<std::vector<double>> h_scalar  (num_GPUs, std::vector<double>(n_bins, 0));


    #pragma omp parallel for num_threads(num_GPUs) // start each GPU running - they each get a single CPU thread as manager
    for (int gpu = 0; gpu < num_GPUs; ++gpu) {
        cudaSetDevice(gpu); // select the current gpu device which will be used for all CUDA calls in this thread (i.e. loop)

        // Calculate the start and end indices for this GPU's workload - there should be no overlap between GPUs and work is equal
        // (technically the final GPU may have slightly less work if total_iterations is not divisible by num_GPUs)
        long long iter_start = gpu * pairs_per_gpu;
        // this is the exclusive end index, so the last GPU may have fewer pairs if total_iterations is not divisible by num_GPUs
        long long iter_end   = std::min(iter_start + pairs_per_gpu, total_iterations); 

        // --- Per-GPU device pointers ---
        double*  d_bins_g;
        double3* d_positions_g;
        double*  d_global_bin_hist_g;
        double3* d_velocities_g  = nullptr;
        double*  d_scalars_g     = nullptr;
        double*  d_global_vsf_g       = nullptr;
        double*  d_global_vsf_parr_g  = nullptr;
        double*  d_global_vsf_perp_g  = nullptr;
        double*  d_global_scalar_g    = nullptr;

        // --- Allocate and copy (same as your existing block, just using _g pointers) ---
        cudaMalloc(&d_positions_g, N * sizeof(double3));
        cudaMemcpy(d_positions_g, positions_double3.data(), N * sizeof(double3), cudaMemcpyHostToDevice);

        cudaMalloc(&d_bins_g, n_bins * sizeof(double));
        cudaMemcpy(d_bins_g, bins.data(), n_bins * sizeof(double), cudaMemcpyHostToDevice);

        cudaMalloc(&d_global_bin_hist_g, n_bins * sizeof(double));
        cudaMemset(d_global_bin_hist_g, 0, n_bins * sizeof(double));

        // Optional allocations based on input files

        //TODO: Change velocities to generic vector
        if (!velocities.empty()) {
            cudaMalloc(&d_velocities_g, N * sizeof(double3));
            cudaMemcpy(d_velocities_g, velocities_double3.data(), N * sizeof(double3), cudaMemcpyHostToDevice);
        }
        if (!scalars.empty()) {
            cudaMalloc(&d_scalars_g, N * sizeof(double));
            cudaMemcpy(d_scalars_g, scalars.data(), N * sizeof(double), cudaMemcpyHostToDevice);
        }

        // if doing the vector structure function, allocate memory for the results on the GPU
        if (do_vsf) {
            cudaMalloc(&d_global_vsf_g, n_bins * sizeof(double));
            cudaMemset(d_global_vsf_g, 0, n_bins * sizeof(double));
        }
        // If projecting vector structure function parallel to the separation vector, allocate memory for the results on the GPU
        if (do_vsf_parr) {
            cudaMalloc(&d_global_vsf_parr_g, n_bins * sizeof(double));
            cudaMemset(d_global_vsf_parr_g, 0, n_bins * sizeof(double));
        }
        // If projecting vector structure function perpendicular to the separation vector, allocate memory for the results on the GPU
        if (do_vsf_perp) {
            cudaMalloc(&d_global_vsf_perp_g, n_bins * sizeof(double));
            cudaMemset(d_global_vsf_perp_g, 0, n_bins * sizeof(double));
        }
        // If doing scalar structure function, allocate memory for the results on the GPU
        if (do_scalar) {
            cudaMalloc(&d_global_scalar_g, n_bins * sizeof(double));
            cudaMemset(d_global_scalar_g, 0, n_bins * sizeof(double));
        }

        // --- Kernel launch config (per GPU) ---
        int blockSize = 256;
        // Get the number of streaming multiprocessors (SMs) on the GPU to determine grid size
        // This is dynamic and will adapt to different GPU architectures (i.e. will work for an old Titan Xp as well as a modern H100)
        int numSM; 
        cudaDeviceGetAttribute(&numSM, cudaDevAttrMultiProcessorCount, 0);
        int gridSize = numSM * 8;
        size_t sharedMemSize = 5 * n_bins * sizeof(double); // (4 histograms + 1 bin_hist)

        std::cout << "GPU " << gpu << ": pairs " << iter_start << " to " << iter_end
                << ", launching " << gridSize << " blocks x " << blockSize << " threads\n";

        vsf_kernel<<<gridSize, blockSize, sharedMemSize>>>(
            N,                   // number of cells
            n_bins,              // number of bins for the structure function
            d_bins_g,            // point to bin edges in global memory
            iter_start,          // start index for this kernel launch
            iter_end,            // end index for this kernel launch
            d_global_bin_hist_g, // point to global histogram for bin counts in global memory
            d_positions_g,       // point to position data in global memory
            d_velocities_g,      // (optional) point to vector data in global memory
            d_scalars_g,         // (optional) point to scalar quantity data in global memory
            d_global_vsf_g,      // (optional) point to global histogram for vector differences in global memory
            d_global_vsf_parr_g, // (optional) point to global histogram for parallel vector differences in global memory
            d_global_vsf_perp_g, // (optional) point to global histogram for perpendicular vector differences in global memory
            d_global_scalar_g    // (optional) point to global histogram for scalar quantity in global memory
        );

        // Check for kernel launch errors
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
            std::cerr << "GPU " << gpu << " kernel launch error: " << cudaGetErrorString(err) << "\n";

        // Wait for GPU to finish before accessing results - this is per CPU thread - each waits for its own GPU
        cudaDeviceSynchronize();

        // Check for kernel execution errors
        err = cudaGetLastError();
        if (err != cudaSuccess)
            std::cerr << "GPU " << gpu << " kernel execution error: " << cudaGetErrorString(err) << "\n";

        // Copy results into per-GPU host buffers
        cudaMemcpy(h_bin_hist[gpu].data(), d_global_bin_hist_g, n_bins * sizeof(double), cudaMemcpyDeviceToHost);
        if (do_vsf)      cudaMemcpy(h_vsf[gpu].data(),      d_global_vsf_g,      n_bins * sizeof(double), cudaMemcpyDeviceToHost);
        if (do_vsf_parr) cudaMemcpy(h_vsf_parr[gpu].data(), d_global_vsf_parr_g, n_bins * sizeof(double), cudaMemcpyDeviceToHost);
        if (do_vsf_perp) cudaMemcpy(h_vsf_perp[gpu].data(), d_global_vsf_perp_g, n_bins * sizeof(double), cudaMemcpyDeviceToHost);
        if (do_scalar)   cudaMemcpy(h_scalar[gpu].data(),   d_global_scalar_g,   n_bins * sizeof(double), cudaMemcpyDeviceToHost);

        // Deallocate GPU memory
        cudaFree(d_positions_g);
        cudaFree(d_bins_g);
        cudaFree(d_global_bin_hist_g);
        if (d_velocities_g)     cudaFree(d_velocities_g);
        if (d_scalars_g)        cudaFree(d_scalars_g);
        if (d_global_vsf_g)     cudaFree(d_global_vsf_g);
        if (d_global_vsf_parr_g) cudaFree(d_global_vsf_parr_g);
        if (d_global_vsf_perp_g) cudaFree(d_global_vsf_perp_g);
        if (d_global_scalar_g)  cudaFree(d_global_scalar_g);

    } // end GPU loop

    // --------------------------------------------
    // CPU reduction: sum per-GPU histograms into final result vectors
    // --------------------------------------------
    std::vector<double> global_bin_hist(n_bins, 0);
    std::vector<double> global_vsf(n_bins, 0);
    std::vector<double> global_vsf_parr(n_bins, 0);
    std::vector<double> global_vsf_perp(n_bins, 0);
    std::vector<double> global_scalar(n_bins, 0);



    // Loop over GPUs and bins to sum the results
    // Note, this will still work if num_GPUs is 1, as it will just sum the single GPU's results into the global vectors
    for (int gpu = 0; gpu < num_GPUs; ++gpu) {
        for (int b = 0; b < n_bins; ++b) {
            global_bin_hist[b] += h_bin_hist[gpu][b];
            if (do_vsf)      global_vsf[b]      += h_vsf[gpu][b];
            if (do_vsf_parr) global_vsf_parr[b] += h_vsf_parr[gpu][b];
            if (do_vsf_perp) global_vsf_perp[b] += h_vsf_perp[gpu][b];
            if (do_scalar)   global_scalar[b]   += h_scalar[gpu][b];
        }
    }


    // --------------------------------------------
    // Normalise VSF by bin counts (if req)
    // --------------------------------------------
    if (norm) {
        std::cout << "Normalising VSF by bin counts...\n";
        for (int i = 0; i < n_bins; ++i) {
            if (global_bin_hist[i] > 0) {
                if (do_vsf) {
                    global_vsf[i] /= global_bin_hist[i]; // normalise by count in each bin
                }
                if (do_vsf_parr) {
                    global_vsf_parr[i] /= global_bin_hist[i]; // normalise by count in each bin
                }
                if (do_vsf_perp) {
                    global_vsf_perp[i] /= global_bin_hist[i]; // normalise by count in each bin
                }
                if (do_scalar) {
                    global_scalar[i] /= global_bin_hist[i]; // normalise by count in each bin
                }
            }
        }
    }

    // --------------------------------------------
    // Save results to binary files
    // --------------------------------------------
    save_to_binary(io_params["output_dir_base"]+io_params["output_dir"]+
                    io_params["output_file_prefix"]+
                    "_"+std::to_string(snap)+"_"+io_params["region"]+
                    "_bin_data.bin", bins);
        
    save_to_binary(io_params["output_dir_base"]+io_params["output_dir"]+
                    io_params["output_file_prefix"]+
                    "_"+std::to_string(snap)+"_"+io_params["region"]+
                    "_counts_data.bin", global_bin_hist);

    if (do_vsf){
        save_to_binary(io_params["output_dir_base"]+io_params["output_dir"]+
                        io_params["output_file_prefix"]+
                        "_"+std::to_string(snap)+"_"+io_params["region"]+
                        "_VSF_data.bin", global_vsf);
    }

    if (do_vsf_parr){
        save_to_binary(io_params["output_dir_base"]+io_params["output_dir"]+
                        io_params["output_file_prefix"]+
                        "_"+std::to_string(snap)+"_"+io_params["region"]+
                        "_VSF_parr_data.bin", global_vsf_parr);
    }

    if (do_vsf_perp){
        save_to_binary(io_params["output_dir_base"]+io_params["output_dir"]+
                        io_params["output_file_prefix"]+
                        "_"+std::to_string(snap)+"_"+io_params["region"]+
                        "_VSF_perp_data.bin", global_vsf_perp);
    }

    if (do_scalar){
        save_to_binary(io_params["output_dir_base"]+io_params["output_dir"]+
                        io_params["output_file_prefix"]+
                        "_"+std::to_string(snap)+"_"+io_params["region"]+
                        "_"+sc_prop+"_data.bin", global_scalar);
    }

    // // --------------------------------------------
    // // Free GPU memory
    // // --------------------------------------------
    // cudaFree(d_positions);
    // cudaFree(d_bins);
    // cudaFree(d_velocities);
    // cudaFree(d_scalars);
    // cudaFree(d_global_bin_hist);
    // cudaFree(d_global_vsf);
    // cudaFree(d_global_vsf_parr);
    // cudaFree(d_global_vsf_perp);
    // cudaFree(d_global_scalar);


    return 0;
}

