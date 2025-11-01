Implementation of Fast Blind Rotation for Bootstrapping FHEs
=====================================

## Fast Blind Rotation for Bootstrapping FHEs
The CHIFHE library contains the implementation of the fully homorphic encryption schemes presented in the paper [Fast Blind Rotation for Bootstrapping FHEs](https://eprint.iacr.org/2023/1564) by using [OpenFHE_v1.1.1](https://github.com/openfheorg/openfhe-development/releases/tag/v1.1.1).

For optimal performance, we also employ an approximate gadget decomposition and provide improved parameter sets as in LMKCDEY (see `binfhecontext.cpp`).
### Requirements
A C++ compiler, the NTL libraries.

## Run the code
1. Configure, build and compile the project.
```
mkdir build
cd build
cmake -DWITH_NTL=ON .. 
make 
```
2. Run the `boolean-xzddf` program in `build/bin/examples/binfhe`.
   
Experimental Result(12th Gen Intel(R) Core(TM) i9-12900H @2.50 GHz and 32 GB RAM, running Ubuntu 20.04.6 LTS):
 ```
100 times of XZDDF      P128G    gate bootstrapping:    2812.5ms
100 times of LMKCDEY    STD128_LMKCDEY   gate bootstrapping:    5000ms
100 times of DM STD128   gate bootstrapping:    6468.75ms
100 times of CGGI       STD128   gate bootstrapping:    6796.87ms
```
We recommend using the following CMake command-line configuration for best performance.
```
cmake -DWITH_NTL=ON  -DNATIVE_SIZE=32 -DWITH_NATIVEOPT=ON -DCMAKE_C_COMPILER=clang-12 -DCMAKE_CXX_COMPILER=clang++-12 -DWITH_OPENMP=OFF -DCMAKE_C_FLAGS="-pthread" -DCMAKE_CXX_FLAGS="-pthread" .. 
```
Experimental Result(12th Gen Intel(R) Core(TM) i9-12900H @2.50 GHz and 32 GB RAM, running Ubuntu 20.04.6 LTS):
```
100 times of XZDDF      P128G    gate bootstrapping:    468.75ms
100 times of LMKCDEY    STD128_LMKCDEY   gate bootstrapping:    578.125ms
100 times of DM STD128   gate bootstrapping:    781.25ms
100 times of CGGI       STD128   gate bootstrapping:    765.625ms
```
