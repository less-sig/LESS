TODO remove

18bb0daf28bfc806649ad406c6638f3b9c07ebb5  PQCsignKAT_105174.rsp
8df56aa1ba520ffac241d8294f12f11ed5accf1d  PQCsignKAT_13940.rsp
c9df47972b63d915e1693ad687722bf46f43e7cd  PQCsignKAT_197315.rsp
8398b442844bfd00f40d329e45a8f8e78fe8cb6e  PQCsignKAT_35074.rsp
5978466ec06c646a027d8305ad232bedc636c2ad  PQCsignKAT_41788.rsp
c3b94175d0367439cba18d6990cf8c1d5a8b8b71  PQCsignKAT_65793.rsp
98019ba6113c9f1743955a0fb120d313838c9cef  PQCsignKAT_97484.rsp

AVX2 
CHES paper:
sign: 117.0
vrfy: 106.8

After Gauß fix
Key generation kCycles (avg,stddev):   725.72,21.88
Signature kCycles (avg,stddev):     134955.52,4338.17
Verification kCycles (avg,stddev):  110472.24,1282.03

After right mul fix 
Key generation kCycles (avg,stddev):   649.47,45.52
kSignature kCycles (avg,stddev):    130633.64,1836.56
Verification kCycles (avg,stddev):  107606.55,1817.26

After normalized_copy_from_generator_non_information_set
Key generation kCycles (avg,stddev):   658.12,42.94
Signature kCycles (avg,stddev):     124807.28,3526.02
Verification kCycles (avg,stddev):  102872.68,4438.76

AVX512 
CHES paper:
sign: 112.2
vrfy:  80.6

After Gauß fix:
Key generation kCycles (avg,stddev):   466.86,19.30
Signature kCycles (avg,stddev):     123132.42,857.23
Verification kCycles (avg,stddev):   85493.60,725.60

After right mul fix 
Key generation kCycles (avg,stddev):   489.47,52.89
Signature kCycles (avg,stddev):     120312.14,5400.90
Verification kCycles (avg,stddev):   82619.34,545.12

The submission for LESS contains the following:

## Reference Implementation

A reference implementation of the signature LESS as a C99 library. The 
reference implementation library does not provide a main() function, and is
intended to be compiled and linked to a binary. The required NIST API is 
present in `Reference_Implementation/include/api.h`. The implementation employ 
a provided generic keccak implementation.

## Optimized Implementation

We provide multiple different optimized versions of the LESS signature scheme.
In total we provide three different optimized implementation: an AVX2, AVX512 
and  NEON. The AVX2 version is optimized for modern Intel CPUs from the Haswell 
generation onwards. The AVX512 implementation targets very modern AMD/Intel 
CPUs from the Zen4/Icelake generation onwards. The NEON version is optimized 
for ARM CPUs like Apple M1. 

The library is organized in the same fashion as the reference one. All 
implementations come with their own optimized SHA3 and SHAKE implementation.

The following libraries/tools are needed:
- `cmake`
- `make`
- `gcc`/`clang`
- `pkg-config`
- `openssl`


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
which translates to:
1 - Enter Utilities/Benchmarking
2 - create a "build" directory and enter it
3 - type cmake ../ to generate the makefiles
4 - type make to compile the codebase

By default, the optimized AVX2 implementation will be compiled. To select the
reference implementation for building, add `-DUSE_REFERENCE=1` to the `cmake`
command in step 3. To enable the ARM Neon implementation you need to add a 
`-DUSE_NEON=1` to the cmake command. The AVX512 implementation is enabled via 
the flag `-DUSE_AVX512=1`. Note it is only possible to pass a single of those 4 
flags (`-DUSE_REFERENCE=1,-DUSE_AVX2=1,-DUSE_AVX512=1,-DUSE_NEON=1`)to cmake.

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

The `KAT_Generation` directory is organized in the same fashion as the
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
-DUSE_AVX512=1
-DUSE_NEON=1
-DUSE_REFERENCE=1
```

to the command in step 3. Only ever add one of the targets. The standard is to 
use the reference implementation. Additionally, you can also enable a `debug` 
build by adding `-DCMAKE_BUILD_TYPE=Debug` to the command in step 3. Furthermore
you can add `-DUSE_SANITIZE=1` to enable memory/pointer sanitation.

A commodity script generating all KATs is provided: `gen_all_kat.sh`. It can be 
run without parameters. KATs will be generated in the KAT directory

NOTE: If you are on a apple system and want to use the memory sanitizer: make 
sure that you use clang, by appending `-DCMAKE_C_COMPILER=clang`. Otherwise 
`LESS` will most likely not build in  debug mode, due to missing support for
the address sanitizer.

Status:
=======
All [proposed improvements](https://eprint.iacr.org/2025/1424.pdf) are implemented.
