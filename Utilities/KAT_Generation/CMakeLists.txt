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

option(USE_SANITIZE "enable compiler memory sanitation" OFF)

# default compilation picks reference codebase
option(USE_REFERENCE "Use reference implementation" ON)

# use AVX2 implementations.
option(USE_AVX2 "Use the AVX2 Optimized Omplementation." OFF)
option(USE_NEON "Use the NEON Optimized Omplementation." OFF)
option(USE_M4 "Use the M4 Optimized Omplementation." OFF)

if(USE_AVX2)
  set(USE_REFERENCE OFF)
endif()

if(USE_NEON)
  set(USE_REFERENCE OFF)
endif()

if(USE_M4)
  set(USE_REFERENCE OFF)
endif()

include(CheckCCompilerFlag)
unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
check_c_compiler_flag(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)

include(CheckIPOSupported)
check_ipo_supported(RESULT lto_supported OUTPUT error)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wpedantic -pedantic -Wuninitialized -Wsign-conversion -Wno-strict-prototypes -Wno-unused-function")

if(USE_SANITZE)
  set(SANITIZE "-fsanitize=address -fsanitize=pointer-compare")
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" UPPER_CMAKE_BUILD_TYPE)
if(UPPER_CMAKE_BUILD_TYPE MATCHES RELEASE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mtune=native -O3 -g0")
  if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
  endif()
  if(lto_supported)
    message(STATUS "IPO / LTO enabled")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
else() # Debug, RelWithDebInfo, etc
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -fsanitize=address")
  message(STATUS "Building in Debug mode!")
endif()

set(CMAKE_C_STANDARD 99)


# add needed optimization flags
# this is generally not needed as `-march=native` is enough to hint the compiler 
# that certain optimizations are available. But not so on nixos.
if(USE_AVX2)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx -mavx2")
endif()

if(USE_NEON)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -flax-vector-conversions")
endif()

message(STATUS "Compilation flags:" ${CMAKE_C_FLAGS})


find_library(KECCAK_LIB keccak)
if(NOT KECCAK_LIB)
 set(STANDALONE_KECCAK 1)
endif()


# selection of specialized compilation units differing between ref and opt
# implementations. 
if(USE_REFERENCE)
  message(STATUS "Using the Reference Implementation")
  set(BASE_DIR ../../Reference_Implementation)
  set(USE_AVX2 OFF)
  set(USE_NEON OFF)
  set(USE_M4 OFF)
elseif(USE_AVX2)
  message(STATUS "Using the AVX2 Implementation")
  set(BASE_DIR ../../Optimized_Implementation/avx2)
  set(USE_NEON OFF)
  set(USE_M4 OFF)
elseif(USE_NEON)
  message(STATUS "Using the NEON Implementation")
  set(BASE_DIR ../../Optimized_Implementation/neon)
  set(USE_AVX2 OFF)
  set(USE_M4 OFF)
elseif(USE_M4)
  message(STATUS "Using the M4 Implementation")
  set(BASE_DIR ../../Optimized_Implementation/m4)
  set(USE_AVX2 OFF)
  set(USE_NEON OFF)
else()
  message( FATAL_ERROR "No implementation specified, aborting." )
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
  set(KECCAK_EXTERNAL_LIB "")
  set(KECCAK_EXTERNAL_ENABLE "")
  list(APPEND COMMON_SOURCES ${BASE_DIR}/lib/keccakf1600.c)
  list(APPEND COMMON_SOURCES ${BASE_DIR}/lib/fips202.c)
else()
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

if(USE_AVX2)
	set(SOURCES ${SOURCES} ${BASE_DIR}/lib/KeccakP-1600-AVX2.s)
	set(HEADERS ${HEADERS} ${BASE_DIR}/include/macro.h)
endif()

if(USE_NEON)
	set(HEADERS ${HEADERS} ${BASE_DIR}/include/macro.h)
endif()

find_package(OpenSSL)

foreach(category RANGE 1 5 2)
    if(category EQUAL 1)
    set(PARAM_TARGETS SHORT_SIG INTERMEDIATE BALANCED)
    else()
    set(PARAM_TARGETS SHORT_SIG BALANCED)
    endif()
    foreach(optimiz_target ${PARAM_TARGETS})
         # settings for benchmarking binary
         set(TARGET_BINARY_NAME LESS_KATgen_cat_${category}_${optimiz_target})
         add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES}
                            ./include/KAT_NIST_rng.h
                            ./lib/KAT_NIST_rng.c
                            ./lib/PQCgenKAT_sign.c)
         target_include_directories(${TARGET_BINARY_NAME} PRIVATE
           ${BASE_DIR}/include
           ./include
           ${OPENSSL_CRYPTO_LIBRARY}
         )
         target_link_libraries(${TARGET_BINARY_NAME} PRIVATE OpenSSL::Crypto ${SANITIZE} ${KECCAK_EXTERNAL_LIB})
         set_target_properties(${TARGET_BINARY_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./bin)
         set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY
            COMPILE_FLAGS "-DCATEGORY_${category} -D${optimiz_target}=1 ${KECCAK_EXTERNAL_ENABLE} ")
    endforeach(optimiz_target)
endforeach(category)
