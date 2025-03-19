cmake_minimum_required(VERSION 3.0)
project(less)
enable_language(C ASM)
set(CMAKE_C_STANDARD 99)
set(CMAKE_VERBOSE_MAKEFILE OFF)

enable_testing()
find_package(OpenSSL)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# needed for nvim
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
if(CMAKE_EXPORT_COMPILE_COMMANDS)
  set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES
          ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()

set(ALLOWED_WARNINGS " -Wno-type-limits")
IF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm")
	set(CMAKE_C_FLAGS_DEBUG "-DDEBUG -Og -ggdb -Wall -Wextra -Wpedantic -fsanitize=address ${ALLOWED_WARNINGS}")
	set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -Wall -Wextra -Wpedantic -march=native -mtune=native -flax-vector-conversions -ftree-vectorize -funroll-loops -flto ${ALLOWED_WARNINGS}")
else()
	set(CMAKE_C_FLAGS_DEBUG "-DDEBUG -O0 -ggdb -Wall -Wextra -Wpedantic -fsanitize=address ${ALLOWED_WARNINGS}")
    set(CMAKE_C_FLAGS_RELEASE "-DNDEBUG -O3 -Wall -Wextra -Wpedantic -mavx2 -mavx -march=native -mtune=native -ftree-vectorize -flto ${ALLOWED_WARNINGS}")
endif()


set(SOURCES
        ${PROJECT_SOURCE_DIR}/lib/codes.c
        ${PROJECT_SOURCE_DIR}/lib/fips202.c
        ${PROJECT_SOURCE_DIR}/lib/keccakf1600.c
        ${PROJECT_SOURCE_DIR}/lib/LESS.c
        ${PROJECT_SOURCE_DIR}/lib/monomial.c
        ${PROJECT_SOURCE_DIR}/lib/rng.c
        ${PROJECT_SOURCE_DIR}/lib/seedtree.c
        ${PROJECT_SOURCE_DIR}/lib/sign.c
        ${PROJECT_SOURCE_DIR}/lib/utils.c
        ${PROJECT_SOURCE_DIR}/lib/sort.c
        ${PROJECT_SOURCE_DIR}/lib/canonical.c
        ${PROJECT_SOURCE_DIR}/lib/transpose.c
)

set(HEADERS
        ${PROJECT_SOURCE_DIR}/include/api.h
        ${PROJECT_SOURCE_DIR}/include/codes.h
        ${PROJECT_SOURCE_DIR}/include/fips202.h
        ${PROJECT_SOURCE_DIR}/include/keccakf1600.h
        ${PROJECT_SOURCE_DIR}/include/sha3.h
        ${PROJECT_SOURCE_DIR}/include/fq_arith.h
        ${PROJECT_SOURCE_DIR}/include/LESS.h
        ${PROJECT_SOURCE_DIR}/include/monomial_mat.h
        ${PROJECT_SOURCE_DIR}/include/parameters.h
        ${PROJECT_SOURCE_DIR}/include/rng.h
        ${PROJECT_SOURCE_DIR}/include/seedtree.h
        ${PROJECT_SOURCE_DIR}/include/utils.h
        ${PROJECT_SOURCE_DIR}/include/sort.h
        ${PROJECT_SOURCE_DIR}/include/transpose.h
)

include_directories(include)

foreach(category "252" "400" "548")
  if(${category} EQUAL "252")
    set(PARAM_TARGETS "192" "68" "45")
  elseif(${category} EQUAL "400")
    set(PARAM_TARGETS "220" "102")
  elseif(${category} EQUAL "548")
    set(PARAM_TARGETS "345" "137")
  endif()

  foreach(optimize_target ${PARAM_TARGETS})
    # settings for benchmarking binary
    set(TARGET_BINARY_NAME LESS_benchmark_cat_${category}_${optimize_target})
    add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES} ${PROJECT_SOURCE_DIR}/lib/bench/less_benchmark.c)
    target_link_libraries(${TARGET_BINARY_NAME} m)
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY=${category} -DTARGET=${optimize_target}")
    target_include_directories(${TARGET_BINARY_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lib/test)
    add_test(${TARGET_BINARY_NAME} ${TARGET_BINARY_NAME})

    # settings for unit tests binary
    set(TARGET_BINARY_NAME LESS_test_cat_${category}_${optimize_target})
    add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES} ${PROJECT_SOURCE_DIR}/lib/test/less_test.c)
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY=${category} -DTARGET=${optimize_target}")
    target_include_directories(${TARGET_BINARY_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lib/test)
    add_test(${TARGET_BINARY_NAME} ${TARGET_BINARY_NAME})

    # KATS generation
    set(TARGET_BINARY_NAME LESS_nist_cat_${category}_${optimize_target})
    add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES}  ${PROJECT_SOURCE_DIR}/lib/nist/KAT_NIST_rng.c ${PROJECT_SOURCE_DIR}/lib/nist/PQCgenKAT_sign.c)
    target_include_directories(${TARGET_BINARY_NAME} PRIVATE ${OPENSSL_INCLUDE_DIR})
    target_link_libraries(${TARGET_BINARY_NAME} PRIVATE OpenSSL::Crypto)
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY=${category} -DTARGET=${optimize_target}")
    add_test(${TARGET_BINARY_NAME} ${TARGET_BINARY_NAME})
  endforeach(optimize_target)
endforeach(category)
