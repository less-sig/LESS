cmake_minimum_required(VERSION 3.0.0)
project(LESS)
enable_language(C ASM)
set(CMAKE_C_STANDARD 99)
set(CMAKE_VERBOSE_MAKEFILE OFF)

# Set a default build type if none was specified
set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

option(USE_SANITIZE "enable compiler memory sanitation" OFF)

# use AVX2 implementations.
option(USE_AVX2 "Use the AVX2 Optimized Omplementation." ON)
option(USE_NEON "Use the NEON Optimized Omplementation." OFF)
# option(USE_M4 "Use the M4 Optimized Omplementation." OFF)

# default compilation picks reference codebase
if(USE_REFERENCE)
  message(STATUS "Using the Reference Implementation")
  set(BASE_DIR ../../Reference_Implementation)
  set(USE_AVX2 OFF)
  set(USE_NEON OFF)
#  set(USE_M4 OFF)
elseif(USE_NEON)
  message(STATUS "Using the NEON Implementation")
  set(BASE_DIR ../../Optimized_Implementation/neon)
  set(USE_AVX2 OFF)
#  set(USE_M4 OFF)
  set(USE_REFERENCE OFF)
#elseif(USE_M4)
#  message(STATUS "Using the M4 Implementation")
#  set(BASE_DIR ../../Optimized_Implementation/m4)
#  set(USE_AVX2 OFF)
#  set(USE_NEON OFF)
#  set(USE_REFERENCE OFF)
elseif(USE_AVX2)
  message(STATUS "Using the AVX2 Implementation")
  set(BASE_DIR ../../Optimized_Implementation/avx2)
  set(USE_NEON OFF)
#  set(USE_M4 OFF)
  set(USE_REFERENCE OFF)
else()
  message( FATAL_ERROR "No implementation specified, aborting." )
endif()

include(CheckCCompilerFlag)
unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
check_c_compiler_flag(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wpedantic -pedantic -Wno-type-limits")

if(USE_SANITIZE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fsanitize=pointer-compare")
endif()

if(CMAKE_BUILD_TYPE MATCHES Release)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mtune=native -O3")
  # Floyd cannot use this as avx512 is not supported by valgring
  #if(COMPILER_SUPPORTS_MARCH_NATIVE)
  #  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
  #endif()
  if(lto_supported)
    message(STATUS "IPO / LTO enabled")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
elseif(CMAKE_BUILD_TYPE MATCHES MinSizeRel)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mtune=native -Os -g0")
  if(lto_supported)
    message(STATUS "IPO / LTO enabled")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
else() # Debug, RelWithDebInfo, etc
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Og ")
  message(STATUS "Building in Debug mode!")
endif()

set(CMAKE_C_STANDARD 99)

if(USE_AVX2)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_AVX2 -mavx -mavx2 -mbmi -mbmi2")
endif()

if(USE_NEON)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_NEON -flax-vector-conversions")
endif()

message(STATUS "Compilation flags:" ${CMAKE_C_FLAGS})


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

set(SOURCES
  ${COMMON_SOURCES}
  ${BASE_DIR}/lib/codes.c
  ${BASE_DIR}/lib/LESS.c
  ${BASE_DIR}/lib/monomial.c
  ${BASE_DIR}/lib/rng.c
  ${BASE_DIR}/lib/seedtree.c
  ${BASE_DIR}/lib/utils.c
  ${BASE_DIR}/lib/sign.c
  ${BASE_DIR}/lib/sort.c
  ${BASE_DIR}/lib/transpose.c
  ${BASE_DIR}/lib/canonical.c
  ${BASE_DIR}/lib/keccakf1600.c
  ${BASE_DIR}/lib/fips202.c
)

if(USE_AVX2)
	set(SOURCES ${SOURCES} ${BASE_DIR}/lib/KeccakP-1600-AVX2.s)
	set(SOURCES ${SOURCES} ${BASE_DIR}/lib/transpose_avx2.c)
	set(HEADERS ${HEADERS} ${BASE_DIR}/include/macro.h)
endif()

if(USE_NEON)
	set(SOURCES ${SOURCES} ${BASE_DIR}/lib/transpose_neon.c)
	set(HEADERS ${HEADERS} ${BASE_DIR}/include/macro.h)
endif()

include_directories(${BASE_DIR}/include)
foreach(category "252" "400" "548")
  if(${category} EQUAL "252")
    set(PARAM_TARGETS "192" "68" "45")
  elseif(${category} EQUAL "400")
    set(PARAM_TARGETS "220" "102")
  elseif(${category} EQUAL "548")
    set(PARAM_TARGETS "345" "137")
  endif()
  
    foreach(optimize_target ${PARAM_TARGETS})
   	set(TARGET_BINARY_NAME LESS_memory_${category}_${optimize_target}_keygen)
    add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES}
                       ./lib/less_keygen.c)
    target_include_directories(${TARGET_BINARY_NAME} PRIVATE
                               ${BASE_DIR}/include
                               ./include)
    target_link_libraries(${TARGET_BINARY_NAME} m ${SANITIZE})
    set_target_properties(${TARGET_BINARY_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./bin)
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY
       COMPILE_FLAGS "-DCATEGORY=${category} -DTARGET=${optimize_target}")

   	set(TARGET_BINARY_NAME LESS_memory_${category}_${optimize_target}_sign)
    add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES} ./lib/less_sign.c)
    target_include_directories(${TARGET_BINARY_NAME} PRIVATE
                               ${BASE_DIR}/include
                               ./include)
    target_link_libraries(${TARGET_BINARY_NAME} m ${SANITIZE})
    set_target_properties(${TARGET_BINARY_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./bin)
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY
       COMPILE_FLAGS "-DCATEGORY=${category} -DTARGET=${optimize_target}")

   	set(TARGET_BINARY_NAME LESS_memory_${category}_${optimize_target}_verify)
    add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES}
                       ./lib/less_verify.c)
    target_include_directories(${TARGET_BINARY_NAME} PRIVATE
                               ${BASE_DIR}/include
                               ./include)
    target_link_libraries(${TARGET_BINARY_NAME} m ${SANITIZE})
    set_target_properties(${TARGET_BINARY_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./bin)
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY
       COMPILE_FLAGS "-DCATEGORY=${category} -DTARGET=${optimize_target}")
  endforeach(optimize_target)
endforeach(category)
