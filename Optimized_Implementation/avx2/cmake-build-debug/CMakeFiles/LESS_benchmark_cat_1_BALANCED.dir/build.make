# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Produce verbose output by default.
VERBOSE = 1

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /nix/store/p4kcxw06la5b9nmaksbss065c7j5dpwc-clion-2024.1.3/clion/bin/cmake/linux/x64/bin/cmake

# The command to remove a file.
RM = /nix/store/p4kcxw06la5b9nmaksbss065c7j5dpwc-clion-2024.1.3/clion/bin/cmake/linux/x64/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/flags.make

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/codes.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/flags.make
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/codes.c.o: /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/codes.c
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/codes.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/codes.c.o"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/codes.c.o -MF CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/codes.c.o.d -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/codes.c.o -c /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/codes.c

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/codes.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/codes.c.i"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/codes.c > CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/codes.c.i

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/codes.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/codes.c.s"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/codes.c -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/codes.c.s

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/fips202.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/flags.make
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/fips202.c.o: /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/fips202.c
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/fips202.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/fips202.c.o"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/fips202.c.o -MF CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/fips202.c.o.d -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/fips202.c.o -c /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/fips202.c

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/fips202.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/fips202.c.i"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/fips202.c > CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/fips202.c.i

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/fips202.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/fips202.c.s"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/fips202.c -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/fips202.c.s

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/keccakf1600.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/flags.make
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/keccakf1600.c.o: /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/keccakf1600.c
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/keccakf1600.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/keccakf1600.c.o"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/keccakf1600.c.o -MF CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/keccakf1600.c.o.d -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/keccakf1600.c.o -c /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/keccakf1600.c

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/keccakf1600.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/keccakf1600.c.i"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/keccakf1600.c > CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/keccakf1600.c.i

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/keccakf1600.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/keccakf1600.c.s"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/keccakf1600.c -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/keccakf1600.c.s

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/LESS.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/flags.make
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/LESS.c.o: /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/LESS.c
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/LESS.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/LESS.c.o"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/LESS.c.o -MF CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/LESS.c.o.d -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/LESS.c.o -c /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/LESS.c

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/LESS.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/LESS.c.i"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/LESS.c > CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/LESS.c.i

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/LESS.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/LESS.c.s"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/LESS.c -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/LESS.c.s

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/monomial.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/flags.make
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/monomial.c.o: /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/monomial.c
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/monomial.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/monomial.c.o"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/monomial.c.o -MF CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/monomial.c.o.d -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/monomial.c.o -c /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/monomial.c

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/monomial.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/monomial.c.i"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/monomial.c > CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/monomial.c.i

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/monomial.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/monomial.c.s"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/monomial.c -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/monomial.c.s

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/rng.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/flags.make
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/rng.c.o: /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/rng.c
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/rng.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/rng.c.o"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/rng.c.o -MF CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/rng.c.o.d -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/rng.c.o -c /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/rng.c

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/rng.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/rng.c.i"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/rng.c > CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/rng.c.i

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/rng.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/rng.c.s"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/rng.c -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/rng.c.s

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/seedtree.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/flags.make
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/seedtree.c.o: /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/seedtree.c
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/seedtree.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/seedtree.c.o"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/seedtree.c.o -MF CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/seedtree.c.o.d -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/seedtree.c.o -c /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/seedtree.c

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/seedtree.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/seedtree.c.i"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/seedtree.c > CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/seedtree.c.i

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/seedtree.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/seedtree.c.s"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/seedtree.c -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/seedtree.c.s

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sign.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/flags.make
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sign.c.o: /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/sign.c
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sign.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building C object CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sign.c.o"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sign.c.o -MF CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sign.c.o.d -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sign.c.o -c /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/sign.c

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sign.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sign.c.i"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/sign.c > CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sign.c.i

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sign.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sign.c.s"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/sign.c -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sign.c.s

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/utils.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/flags.make
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/utils.c.o: /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/utils.c
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/utils.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/utils.c.o"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/utils.c.o -MF CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/utils.c.o.d -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/utils.c.o -c /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/utils.c

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/utils.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/utils.c.i"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/utils.c > CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/utils.c.i

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/utils.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/utils.c.s"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/utils.c -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/utils.c.s

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/canonical.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/flags.make
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/canonical.c.o: /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/canonical.c
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/canonical.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building C object CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/canonical.c.o"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/canonical.c.o -MF CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/canonical.c.o.d -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/canonical.c.o -c /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/canonical.c

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/canonical.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/canonical.c.i"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/canonical.c > CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/canonical.c.i

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/canonical.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/canonical.c.s"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/canonical.c -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/canonical.c.s

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sort.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/flags.make
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sort.c.o: /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/sort.c
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sort.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building C object CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sort.c.o"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sort.c.o -MF CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sort.c.o.d -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sort.c.o -c /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/sort.c

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sort.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sort.c.i"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/sort.c > CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sort.c.i

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sort.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sort.c.s"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/sort.c -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sort.c.s

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/KeccakP-1600-AVX2.s.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/flags.make
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/KeccakP-1600-AVX2.s.o: /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/KeccakP-1600-AVX2.s
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building ASM object CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/KeccakP-1600-AVX2.s.o"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(ASM_DEFINES) $(ASM_INCLUDES) $(ASM_FLAGS) -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/KeccakP-1600-AVX2.s.o -c /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/KeccakP-1600-AVX2.s

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/KeccakP-1600-AVX2.s.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing ASM source to CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/KeccakP-1600-AVX2.s.i"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(ASM_DEFINES) $(ASM_INCLUDES) $(ASM_FLAGS) -E /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/KeccakP-1600-AVX2.s > CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/KeccakP-1600-AVX2.s.i

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/KeccakP-1600-AVX2.s.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling ASM source to assembly CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/KeccakP-1600-AVX2.s.s"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(ASM_DEFINES) $(ASM_INCLUDES) $(ASM_FLAGS) -S /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/KeccakP-1600-AVX2.s -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/KeccakP-1600-AVX2.s.s

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/bench/less_benchmark.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/flags.make
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/bench/less_benchmark.c.o: /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/bench/less_benchmark.c
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/bench/less_benchmark.c.o: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building C object CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/bench/less_benchmark.c.o"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/bench/less_benchmark.c.o -MF CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/bench/less_benchmark.c.o.d -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/bench/less_benchmark.c.o -c /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/bench/less_benchmark.c

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/bench/less_benchmark.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/bench/less_benchmark.c.i"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/bench/less_benchmark.c > CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/bench/less_benchmark.c.i

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/bench/less_benchmark.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/bench/less_benchmark.c.s"
	/nix/store/wk45q17zzwwpsgary6hfc368567kmwiy-clang-wrapper-17.0.6/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/lib/bench/less_benchmark.c -o CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/bench/less_benchmark.c.s

# Object files for target LESS_benchmark_cat_1_BALANCED
LESS_benchmark_cat_1_BALANCED_OBJECTS = \
"CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/codes.c.o" \
"CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/fips202.c.o" \
"CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/keccakf1600.c.o" \
"CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/LESS.c.o" \
"CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/monomial.c.o" \
"CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/rng.c.o" \
"CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/seedtree.c.o" \
"CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sign.c.o" \
"CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/utils.c.o" \
"CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/canonical.c.o" \
"CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sort.c.o" \
"CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/KeccakP-1600-AVX2.s.o" \
"CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/bench/less_benchmark.c.o"

# External object files for target LESS_benchmark_cat_1_BALANCED
LESS_benchmark_cat_1_BALANCED_EXTERNAL_OBJECTS =

LESS_benchmark_cat_1_BALANCED: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/codes.c.o
LESS_benchmark_cat_1_BALANCED: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/fips202.c.o
LESS_benchmark_cat_1_BALANCED: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/keccakf1600.c.o
LESS_benchmark_cat_1_BALANCED: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/LESS.c.o
LESS_benchmark_cat_1_BALANCED: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/monomial.c.o
LESS_benchmark_cat_1_BALANCED: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/rng.c.o
LESS_benchmark_cat_1_BALANCED: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/seedtree.c.o
LESS_benchmark_cat_1_BALANCED: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sign.c.o
LESS_benchmark_cat_1_BALANCED: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/utils.c.o
LESS_benchmark_cat_1_BALANCED: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/canonical.c.o
LESS_benchmark_cat_1_BALANCED: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/sort.c.o
LESS_benchmark_cat_1_BALANCED: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/KeccakP-1600-AVX2.s.o
LESS_benchmark_cat_1_BALANCED: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/lib/bench/less_benchmark.c.o
LESS_benchmark_cat_1_BALANCED: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/build.make
LESS_benchmark_cat_1_BALANCED: CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Linking C executable LESS_benchmark_cat_1_BALANCED"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/build: LESS_benchmark_cat_1_BALANCED
.PHONY : CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/build

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/cmake_clean.cmake
.PHONY : CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/clean

CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/depend:
	cd /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2 /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2 /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug /home/duda/Downloads/crypto/schemes/LESS/Optimized_Implementation/avx2/cmake-build-debug/CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/LESS_benchmark_cat_1_BALANCED.dir/depend

