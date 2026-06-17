NVCC = nvcc
CXX = g++

NVCCFLAGS = -O3 \
    -gencode arch=compute_61,code=sm_61 \
    -gencode arch=compute_89,code=sm_89 \
    -gencode arch=compute_90,code=sm_90 \
    -gencode arch=compute_90,code=compute_90 \
    -Xcompiler -fopenmp

CXXFLAGS = -O3 -fopenmp

LDFLAGS = -lgomp

GPU_SRCS = struc_func.cu struc_func_functions.cu
CPU_SRCS = struc_func.cpp struc_func_functions.cpp

GPU_TARGET = sf_cuda.out
CPU_TARGET = sf_cpu.out

# Default target - builds both
all: $(GPU_TARGET) $(CPU_TARGET)

$(GPU_TARGET): $(GPU_SRCS)
	$(NVCC) $(NVCCFLAGS) $(GPU_SRCS) -o $(GPU_TARGET) $(LDFLAGS)

$(CPU_TARGET): $(CPU_SRCS)
	$(CXX) $(CXXFLAGS) $(CPU_SRCS) -o $(CPU_TARGET) $(LDFLAGS)

clean:
	rm -f $(GPU_TARGET) $(CPU_TARGET)