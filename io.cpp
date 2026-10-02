struct BinningConfig {
    int n_bins;
    double bin_min;
    double bin_max;
    bool logarithmic;
};

struct RuntimeConfig {
    std::string backend;
    int cpu_threads;
    int gpu_count;
    uint64_t random_seed;
};

struct AnalysisConfig {
    BinningConfig binning;
    bool normalise;
    std::string method;
};
