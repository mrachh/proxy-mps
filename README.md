# proxy-mps
This repository contains codes for doing multiple particle scattering using proxy-surfaces

## Installing this package

To compile from git, clone repository with the submodules

    git clone --recurse-submodules https://github.com/mrachh/proxy-mps.git

and run startup.m in the install directory.

This will download chunkIE (including its submodules FLAM and fmm2d), and fmm3dbie 
(including its submodule fmm3dbie), include chunkIE, FLAM in the MATLAB path, 
and generate the relevant mex files for fmm2d, FMM3D, and fmm3dbie mex files
if a fortran compiler exists.

The layered medium codes require an additional dependency which will be added soon.

