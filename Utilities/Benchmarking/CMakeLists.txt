cmake_minimum_required(VERSION 3.9.4)

project(LESS C)

# Set a default build type if none was specified
set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wpedantic -pedantic -Wuninitialized -Wsign-conversion -Wno-strict-prototypes")

include(CheckCCompilerFlag)
unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
check_c_compiler_flag(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)

include(CheckIPOSupported)
check_ipo_supported(RESULT lto_supported OUTPUT error)

set(SANITIZE "")
# set(SANITIZE "-fsanitize=address -fsanitize=pointer-compare")

if(CMAKE_BUILD_TYPE MATCHES Release)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mtune=native -O3 -g0")
  if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
  endif()
  if(lto_supported)
    message(STATUS "IPO / LTO enabled")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
elseif(CMAKE_BUILD_TYPE MATCHES MinSizeRel)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mtune=native -Os -g0")
  if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
  endif()
  if(lto_supported)
    message(STATUS "IPO / LTO enabled")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
else() # Debug, RelWithDebInfo, etc
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fsanitize=pointer-compare")
  message(STATUS "Building in Debug mode!")
endif()

# TODO: NIST seems to mandate "ANSI C" for submissions.
# GCC's "-ansi" flag is equivalent to -std=c00, but current code is at least C99
set(CMAKE_C_STANDARD 99)

message(STATUS "Compilation flags:" ${CMAKE_C_FLAGS})

option(COMPRESS_CMT_COLUMNS "Enable COMPRESS_CMT_COLUMNS to compress commitment in SG and VY before hashing (reduces SHA-3 permutaitons by 12.5%, but adds overhead of encoding)" OFF)
if(COMPRESS_CMT_COLUMNS)
message(STATUS "COMPRESS_CMT_COLUMNS is enabled")
add_definitions(-DCOMPRESS_CMT_COLUMNS)
else()
message(STATUS "COMPRESS_CMT_COLUMNS is disabled")
endif()
unset(COMPRESS_CMT_COLUMNS CACHE)


find_library(KECCAK_LIB keccak)
if(NOT KECCAK_LIB)
 set(STANDALONE_KECCAK 1)
endif()

# selection of specialized compilation units differing between ref and opt implementations.
option(AVX2_OPTIMIZED "Use the AVX2 Optimized Omplementation. If not set the Reference Implementation will be used." OFF)

# default compilation picks reference codebase
if(AVX2_OPTIMIZED)
  message(STATUS "Using the AVX2 Optimized Implementation")
  set(BASE_DIR ../../Additional_Implementations/AVX2)
else()
  message(STATUS "Using the Reference Implementation")
  set(BASE_DIR ../../Reference_Implementation)
endif()

set(HEADERS
  ${BASE_DIR}/include/api.h
  ${BASE_DIR}/include/codes.h
  ${BASE_DIR}/include/fips202.h
  ${BASE_DIR}/include/fq_arith.h
  ${BASE_DIR}/include/keccakf1600.h
  ${BASE_DIR}/include/LESS.h
  ${BASE_DIR}/include/monomial_mat.h
  ${BASE_DIR}/include/parameters.h
  ${BASE_DIR}/include/rng.h
  ${BASE_DIR}/include/seedtree.h
  ${BASE_DIR}/include/sha3.h
  ${BASE_DIR}/include/utils.h
)

if(STANDALONE_KECCAK)
  message(STATUS "Employing standalone SHA-3")
  set(KECCAK_EXTERNAL_LIB "")
  set(KECCAK_EXTERNAL_ENABLE "")
  list(APPEND COMMON_SOURCES ${BASE_DIR}/lib/keccakf1600.c)
  list(APPEND COMMON_SOURCES ${BASE_DIR}/lib/fips202.c)
else()
  message(STATUS "Employing libkeccak")
  set(KECCAK_EXTERNAL_LIB keccak)
  set(KECCAK_EXTERNAL_ENABLE "-DSHA_3_LIBKECCAK")
endif()


set(SOURCES
  ${COMMON_SOURCES}
  ${BASE_DIR}/lib/codes.c
  ${BASE_DIR}/lib/LESS.c
  ${BASE_DIR}/lib/monomial.c
  ${BASE_DIR}/lib/rng.c
  ${BASE_DIR}/lib/seedtree.c
  ${BASE_DIR}/lib/utils.c
  ${BASE_DIR}/lib/sign.c
)

if(AVX2_OPTIMIZED)
  set(SOURCES ${SOURCES} ${BASE_DIR}/lib/avx2_table.c)
  set(HEADERS ${HEADERS} ${BASE_DIR}/include/avx2_macro.h)
endif()

foreach(category RANGE 1 5 2)
  if(category EQUAL 1)
  set(PARAM_TARGETS SHORT_SIG INTERMEDIATE BALANCED)
  else()
  set(PARAM_TARGETS SHORT_SIG BALANCED)
  endif()
  foreach(optimiz_target ${PARAM_TARGETS})
    # settings for benchmarking binary
    set(TARGET_BINARY_NAME LESS_benchmark_cat_${category}_${optimiz_target})
    add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES}
                       ./lib/less_benchmark.c)
    target_include_directories(${TARGET_BINARY_NAME} PRIVATE
                               ${BASE_DIR}/include
                               ./include)
    target_link_libraries(${TARGET_BINARY_NAME} m ${SANITIZE} ${KECCAK_EXTERNAL_LIB})
    set_target_properties(${TARGET_BINARY_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./bin)
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY
       COMPILE_FLAGS "-DCATEGORY_${category}=1 -D${optimiz_target}=1 ${KECCAK_EXTERNAL_ENABLE} ")
    # settings for unit tests binary
    set(TARGET_BINARY_NAME LESS_test_cat_${category}_${optimiz_target})
    add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES}
                       ./lib/less_test.c)
    target_include_directories(${TARGET_BINARY_NAME} PRIVATE
                               ${BASE_DIR}/include
                               ./include)
    target_link_libraries(${TARGET_BINARY_NAME} m ${SANITIZE} ${KECCAK_EXTERNAL_LIB})
    set_target_properties(${TARGET_BINARY_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./bin)
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY
       COMPILE_FLAGS "-DCATEGORY_${category} -D${optimiz_target}=1 ${KECCAK_EXTERNAL_ENABLE} ")
  endforeach(optimiz_target)
endforeach(category)
