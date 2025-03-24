The submission for LESS contains the following:

• Reference_Implementation

A reference implementation of the signature LESS as a C99 library. The 
reference implementation library does not provide a main() function, and is
intended to be compiled and linked to a binary. The required NIST API is 
present in `Reference_Implementation/include/api.h`. The implementation employ 
a provided generic keccak implementation.

• Optimized_Implementation

We provide multiple different optimized versions of the LESS signature scheme.
In total we provide two different optimized implementation: an AVX2, NEON.
The AVX2 version is optimized for modern Intel CPUs from the Haswell generation.
The NEON version is optimized for ARM CPUs like Apple M1. 

The library is organized in the same fashion as the reference one, except there 
is no support for using libkeccak. All implementations come with their own 
optimized SHA3 and SHAKE implementation.

The following libraries/tools are needed:
- `cmake`
- `make`
- `gcc`/`clang`
- `openssl`
- `pkg-config`


Utilities
===========

This directory provides CMake based building facilities to generate executable 
files either performing a benchmark of LESS, or generating the Known Answer 
Test files in the KAT folder, as well as a script used for parameter generation.

Benchmarking:
=============

To build the benchmarking binaries, run the following:
```bash
cd Utilities/Benchmarking
mkdir build && cd build
cmake .. 
make 
```
1 - Enter Utilities/Benchmarking
2 - create a "build" directory and enter it
3 - type cmake ../ to generate the makefiles
4 - type make to compile the codebase

By default, the optimized AVX2 implementation will be compiled. To select the
reference implementation for building, add `-DUSE_REFERENCE=1` to the `cmake`
command in step 3. To enable the ARM neon implementation you need to add a 
`-DUSE_NEON=1`.

Once build you can either run each binary like:
```bash
./LESS_benchmark_cat_252_45
```
or you can use the `bench.sh` script, which will automatically generate a markdown
table with the results. Additionally you can pass the flag `-DUSE_SANITIZER=1` to 
the cmake command to check for memory errors. Note, this slows the program down.

NOTE: If you run the benchmarks on an ARM based Apple computer, make sure that
    you run the binaries/scripts with root rights.

KAT Generation:
--------------

The KAT_Generation directory is organized in the same fashion as the
Benchmarking one. The same compilation procedure enacted for the benchmarking 
binary will generate all the executable files needed to re-generate the Known 
Answer Tests files. Specifically,
```bash
cd Utilities/KAT_Generation
mkdir build && cd build
cmake .. 
make 
```

By default, the reference implementation will be compiled. To select the 
optimized implementation for building, add either 
```
-DUSE_AVX2=1
-DUSE_NEON=1
-DUSE_REFERENCE=1
```

to the command in step 3. Only ever add one of the targets. The standard is to 
use the reference implementation. Additionally, you can also enable a `debug` 
build by adding `-DCMAKE_BUILD_TYPE=Debug` to the command in step 3. Furthermore
you can add `-DUSE_SANITIZE=1` to enable memory/pointer sanitation.

A commodity script generating all KATs is provided: `gen_all_kat.sh`. It can be 
run without parameters. KATs will be generated in the KAT directory

NOTE: If you are on a apple system: make sure that you use clang, by appending
    `-DCMAKE_C_COMPILER=clang`. Otherwise `LESS` will most likely not build in 
    debug mode, due to missing support for the address sanitizer.


# Runtime Results:

Ryzen7600X:

AVX512 (18.03 avx512):
build/LESS_benchmark_cat_252_192 | 0.86582 | 166.123 | 134.578
build/LESS_benchmark_cat_252_45  | 4.44392 | 37.4373 | 35.77
build/LESS_benchmark_cat_252_68  | 2.01583 |  58.278 | 53.103
build/LESS_benchmark_cat_400_102 | 5.92569 | 203.962 | 183.333
build/LESS_benchmark_cat_400_220 | 2.24789 | 442.648 | 372.883
build/LESS_benchmark_cat_548_137 | 12.4027 | 573.668 | 529.82
build/LESS_benchmark_cat_548_345 |  4.7713 | 1456.4  | 1218.53

AVX2 (18.03 avx512):
build/LESS_benchmark_cat_252_192 | 0.76152 | 178.855 | 169.366
build/LESS_benchmark_cat_252_45  | 3.64251 | 40.1315 | 45.6199
build/LESS_benchmark_cat_252_68  | 1.70032 | 62.0678 | 67.1031
build/LESS_benchmark_cat_400_102 | 4.977   | 222.905 | 227.451
build/LESS_benchmark_cat_400_220 | 1.94111 | 485.961 | 458.293
build/LESS_benchmark_cat_548_137 | 10.7414 | 607.812 | 628.359
build/LESS_benchmark_cat_548_345 | 4.18581 | 1542.57 | 1525.19


AVX2 (18.03 main):
build/LESS_benchmark_cat_252_192 | 0.86415 | 198.058 | 186.588
build/LESS_benchmark_cat_252_45  | 4.44641 | 45.3128 | 48.4815
build/LESS_benchmark_cat_252_68  | 2.01227 | 69.6692 | 71.9864
build/LESS_benchmark_cat_400_102 | 6.07455 | 254.833 | 252.697
build/LESS_benchmark_cat_400_220 | 2.52635 | 548.314 | 514.209
build/LESS_benchmark_cat_548_137 | 12.5745 | 702.384 | 718.763
build/LESS_benchmark_cat_548_345 | 4.69773 | 1751.78 | 1740.79

