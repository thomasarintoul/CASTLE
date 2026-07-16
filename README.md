# `CASTLE`: Characterisation of Astrophysical Structure and Turbulence in Localised Environments
Written by: Thomas Rintoul - PhD Student in Computational Astrophysics at Cardiff University.

## What is CASTLE?

CASTLE is designed to calculate structure functions typically used in fields such as astrophysics. I initially developed it to calculate velocity structure functions for cosmological zoom-in simulations of massive galaxies in my work on the Circumgalactic Medium. It is however, generalisable to any 1D, 2D or 3D dataset provided the elements considered are in some sort of Cartesian co-ordinate system.
It has been extended to allow "scalar structure functions" for properties such as Temperature to be calculated (see section below on the underlying mathematics for details).

This code is written primarily in modern C++, with the option to use CUDA enabled (i.e. Nvidia) GPU(s) to speed up computation.
The core structure function calculation can be run on a single CPU core, multiple CPUs, a single GPU or multiple GPUs.

This has been run on up to 20 CPU cores without issue and with reasonable scaling (test data to be added in future)
For CUDA runs, the code has been tested up to 8 GPUs on a single node. I cannot guarantee multi-node support.

## Mathematics of Structure Functions
(to be added)

## Code performance
(to be added)

## Planned further developments (in no particular order)
- New method for passing analysis and io parameters
- Dynamic handling of array size to reduce required data preparation for 2D data
- Allow structure functions for an arbitrary number of vector and scalar quantities simultaneously to reduce repetition of binning calculation
- Python front-end
- Support for non-Nvidia GPUs
- Multi-node support for HPC use
- Unit tests

## Citation Policy
If you use this code in any published research work, please link to this repository pending a full public release, at which point a DOI will be generated.

## Contributions
I welcome contributions to this project, particularly if you wish to assist expanding the code to work on AMD/Intel GPUs as I do not have ready access to those for testing. Contact me via email: rintoulta {at} cardiff dot ac dot uk
