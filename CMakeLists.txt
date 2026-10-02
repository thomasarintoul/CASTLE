cmake_minimum_required(VERSION 3.18)
project(structure_functions LANGUAGES CXX)

find_package(OpenMP REQUIRED COMPONENTS CXX)

add_executable(sf_cpu struc_func.cpp struc_func_functions.cpp)
target_link_libraries(sf_cpu PRIVATE OpenMP::OpenMP_CXX)

include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
  enable_language(CUDA)
  add_executable(sf_cuda struc_func.cu struc_func_functions.cu)
  set_target_properties(sf_cuda PROPERTIES CUDA_ARCHITECTURES "61;89;90")
  target_link_libraries(sf_cuda PRIVATE OpenMP::OpenMP_CXX)
  target_compile_options(sf_cuda PRIVATE
      $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-fopenmp>)
endif()