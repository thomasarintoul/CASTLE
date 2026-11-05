#pragma once
#ifdef __CUDACC__
  #define HD __host__ __device__
  #define H  __host__
#else
  #define HD
  #define H
#endif