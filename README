The submission for LESS contains the following:

• Reference_Implementation

A reference implementation of the signature LESS as a C11 library. The 
reference implementation library does not provide a main() function, and is
intended to be compiled and linked to a binary. The required NIST API is 
present in Reference_Implementation/include/api.h,

The implementation may either use libkeccak, if installed and available in 
system-wide searched paths (e.g., `/usr/local/include` for the headerfiles and
`/usr/local/lib` for the static library), or employ the provided fallback 
implementation for systems where libkeccak is not  available. 

If the library is not compiled via the supplied CMake flow, defining 
SHA_3_LIBKECCAK during compilation enables compilation of the codebase against 
an available libkeccak installation.

• Optimized_Implementation

We provide multiple different optimized versions of the LESS signature scheme.
In total we provide two different optimized implementation: an AVX2, NEON.
The AVX2 version is optimized for modern Intel CPUs from the Haswell generation.
The NEON version is optimized for ARM CPUs like Apple M1. 

The library is organized in the same fashion as the reference one, except there 
is no support for using libkeccak. All implementations come with their own 
optimized SHA3 and SHAKE implementation.


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
`-DUSE_NEON`.

Once build you can either run each binary like:
```bash
./bin/LESS_benchmark_cat_1_BALANCED
```
or you can use the `bench.sh` script, which will automatically generate a markdown
table with the results.

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
    `-DCMAKE_C_COMPILER=clang`. Otherwise `LESS` will likely not build in 
    debug mode, due to missing support for the address sanitizer.


Runtime Results:
===============

AVX2: (Ryzen 7600X)
|Name            | KeyGen | Sign | Vrfy|
|----------------|--------|------|-----|
|cat_1_BALANCED | 2040.94 | 536298.51 | 555500.31|
|cat_1_INTERMEDIATE | 5629.56 | 543080.11 | 561288.96|
|cat_1_SHORT_SIG | 12471.49 | 433841.97 | 447269.21|
|cat_3_BALANCED | 7060.29 | 5112977.70 | 5333601.13|
|cat_3_SHORT_SIG | 12531.23 | 6218630.60 | 6526643.77|
|cat_5_BALANCED | 14841.81 | 19773892.52 | 20679633.0|
|cat_5_SHORT_SIG | 26316.45 | 13868415.53 | 14488839.30|

REF: (Ryzen 7600X)
|Name            | KeyGen | Sign | Vrfy|
|----------------|--------|------|-----|
|cat_1_BALANCED | 4340.04 | 993920.92 | 984969.12|
|cat_1_INTERMEDIATE | 11230.76 | 971417.44 | 975866.01|
|cat_1_SHORT_SIG | 26883.85 | 815910.79 | 829120.37|
|cat_3_BALANCED | 13305.64 | 10413855.97 | 10598378.74|
|cat_3_SHORT_SIG | 26165.98 | 12535153.75 | 12616438.38|
|cat_5_BALANCED | 35625.37 | 50556007.65 | 51115479.17|
|cat_5_SHORT_SIG | 69411.23 | 33139406.09 | 33322599.90|

NEON: (M1 MacBook Pro)
|Name            | KeyGen | Sign | Vrfy|
|----------------|--------|------|-----|
|cat_1_BALANCED | 1889.76 | 539753.72 | 557659.60|
|cat_1_INTERMEDIATE | 5214.91 | 532322.96 | 551006.21|
|cat_1_SHORT_SIG | 11877.07 | 431441.03 | 446850.50|
|cat_3_BALANCED | 6135.51 | 5264174.06 | 5393275.09|
|cat_3_SHORT_SIG | 11612.32 | 6206678.59 | 6369067.37|
|cat_5_BALANCED | 12618.70 | 19706940.62 | 20193936.79|
|cat_5_SHORT_SIG | 24756.96 | 13984785.27 | 13565152.10|
