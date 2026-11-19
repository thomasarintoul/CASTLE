# Define minimum CMake version and project description
cmake_minimum_required(VERSION 3.18)
project(CASTLE VERSION 0.1.0 LANGUAGES CXX)

# Build options
option(BUILD_CPU "Build CPU implementation" ON)
option(BUILD_CUDA "Build CUDA implementation" OFF)
option(BUILD_HIP "Build HIP implementation" OFF)

# Find OpenMP for parallelisation with C++
find_package(OpenMP REQUIRED COMPONENTS CXX)

### CPU Implementation
if(BUILD_CPU)
  message(STATUS "Building CPU implementation")

  # Define CPU Executable and require C++17 standard
  add_executable(castle_cpu src/main.cpp)
  target_compile_features(castle_cpu PRIVATE cxx_std_17)
  target_link_libraries(castle_cpu PRIVATE OpenMP::OpenMP_CXX)
endif()


if(BUILD_CUDA)
  message(STATUS "Building CUDA implementation")

  # CUDA Implementation
  include(CheckLanguage)
  check_language(CUDA)


  # Default to the host GPU architecture, unless overridden on the CLI
  ### Can override with -DCMAKE_CUDA_ARCHITECTURES=xx on command line
  ### Examples:
  ### 61 for Nvidia Titan Xp
  ### 89 for Nvidia L40S
  ### 90 for Nvidia H100
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES "native" CACHE STRING
        "CUDA architectures to build for")
  endif()

  # Assuming that a CUDA compiler is available
  if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)

    # Isolate CUDA source properties so castle_cpu remains pure C++
    configure_file(src/main.cpp "${CMAKE_CURRENT_BINARY_DIR}/main_cuda.cpp" COPYONLY)
    set_source_files_properties("${CMAKE_CURRENT_BINARY_DIR}/main_cuda.cpp" PROPERTIES LANGUAGE CUDA)

    add_executable(castle_cuda "${CMAKE_CURRENT_BINARY_DIR}/main_cuda.cpp")

    set_target_properties(castle_cuda PROPERTIES
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED YES)

    target_include_directories(castle_cuda PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")

    target_compile_definitions(castle_cuda PRIVATE
        TOML_RETURN_BOOL_FROM_FOR_EACH_BROKEN=1
        TOML_RETURN_BOOL_FROM_FOR_EACH_BROKEN_ACKNOWLEDGED)

    target_compile_options(castle_cuda PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)

    target_link_libraries(castle_cuda PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif()


### HIP Implementation
if(BUILD_HIP)
  message(STATUS "Building HIP implementation")
  ## Complete when HIP port is written
  message(WARNING "HIP implementation is not yet implemented. Please set BUILD_HIP=OFF.")
endif()

### Tests
## Complete when tests are written

### Benchmarks
## Complete when benchmarks are written

## Note on usability - could make a CMakePresets.json file to make it easier to build different versions of CASTLE