cmake_minimum_required(VERSION 3.22)
project(less)
enable_language(C ASM)
set(CMAKE_C_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(OpenSSL)

# needed for nvim
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
if(CMAKE_EXPORT_COMPILE_COMMANDS)
  set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES
          ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()

set(ALLOWED_WARNINGS " -Wno-gnu-binary-literal -Wno-unused-function -Wno-type-limits")

IF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm") 
	set(CMAKE_C_FLAGS_DEBUG "-g -DDEBUG -O0  -march=native -flax-vector-conversions ${ALLOWED_WARNINGS}")
	set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -Wall -Wextra -Wpedantic -march=native -mtune=native -flax-vector-conversions -ftree-vectorize -funroll-loops -fomit-frame-pointer -fno-stack-protector -flto ${ALLOWED_WARNINGS}")
else()
	set(CMAKE_C_FLAGS_DEBUG "-DDEBUG -O0 -ggdb -Wall -Wpedantic -Wextra -fsanitize=address ${ALLOWED_WARNINGS}")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -Wall -Wextra -Wpedantic -mavx2 -mavx -march=native -ftree-vectorize -fomit-frame-pointer -flto ${ALLOWED_WARNINGS}")
endif()


set(SOURCES
        ${PROJECT_SOURCE_DIR}/lib/codes.c
        ${PROJECT_SOURCE_DIR}/lib/fips202.c
        ${PROJECT_SOURCE_DIR}/lib/keccakf1600.c
        ${PROJECT_SOURCE_DIR}/lib/LESS.c
        ${PROJECT_SOURCE_DIR}/lib/LESS_without_tree.c
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

foreach(category RANGE 1 5 2)
  if(category EQUAL 1)
    set(PARAM_TARGETS SHORT_SIG INTERMEDIATE BALANCED)
  else()
    set(PARAM_TARGETS SHORT_SIG BALANCED)
  endif()
  foreach(optimize_target ${PARAM_TARGETS})
    # settings for benchmarking binary
    set(TARGET_BINARY_NAME LESS_benchmark_cat_${category}_${optimize_target})
    add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES} ${PROJECT_SOURCE_DIR}/lib/bench/less_benchmark.c)
    target_link_libraries(${TARGET_BINARY_NAME} m)
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_${category}=1 -D${optimize_target}=1 ${KECCAK_EXTERNAL_ENABLE} ")

    # settings for unit tests binary
    set(TARGET_BINARY_NAME LESS_test_cat_${category}_${optimize_target})
    add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES} ${PROJECT_SOURCE_DIR}/lib/test/less_test.c)
    target_link_libraries(${TARGET_BINARY_NAME} m)
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_${category} -D${optimize_target}=1 ${KECCAK_EXTERNAL_ENABLE} ")

    # KATS generation
    set(TARGET_BINARY_NAME LESS_nist_cat_${category}_${optimize_target})
    add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES}  ${PROJECT_SOURCE_DIR}/lib/nist/KAT_NIST_rng.c ${PROJECT_SOURCE_DIR}/lib/nist/PQCgenKAT_sign.c)
    target_link_libraries(${TARGET_BINARY_NAME} PRIVATE OpenSSL::Crypto)
    include_directories(${OPENSSL_INCLUDE_DIR})
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_${category} -D${optimize_target}=1 ${KECCAK_EXTERNAL_ENABLE} ")
  endforeach(optimize_target)
endforeach(category)

# canonical forms test
add_executable(test_canonical_form ${PROJECT_SOURCE_DIR}/lib/test/test_canonical.c ${SOURCES})
target_link_libraries(test_canonical_form PRIVATE OpenSSL::Crypto)
set_property(TARGET test_canonical_form APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_1 -DBALANCED=1")

add_executable(inv_table ${PROJECT_SOURCE_DIR}/lib/test/inv_table.c)
target_link_libraries(inv_table PRIVATE OpenSSL::Crypto)
set_property(TARGET inv_table APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_1 -DBALANCED=1")

add_executable(arith ${PROJECT_SOURCE_DIR}/lib/test/arith.c)
target_link_libraries(arith PRIVATE OpenSSL::Crypto)
set_property(TARGET arith APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_1 -DBALANCED=1")

add_executable(test_sort ${PROJECT_SOURCE_DIR}/lib/test/sort.c ${SOURCES})
target_link_libraries(test_sort PRIVATE OpenSSL::Crypto)
set_property(TARGET test_sort APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_1 -DBALANCED=1")
# set_property(TARGET test_sort APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_1 -DBALANCED=1")

# canonical forms bench
add_executable(bench_cf ${PROJECT_SOURCE_DIR}/lib/bench/bench_cf.c ${SOURCES})
target_link_libraries(bench_cf PRIVATE OpenSSL::Crypto)
set_property(TARGET bench_cf APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_1 -DBALANCED=1")

add_executable(bench_rref ${PROJECT_SOURCE_DIR}/lib/bench/bench_rref.c ${SOURCES})
target_link_libraries(bench_rref PRIVATE OpenSSL::Crypto)
set_property(TARGET bench_rref APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_1 -DBALANCED=1")

add_executable(bench_sorting ${PROJECT_SOURCE_DIR}/lib/bench/bench_sorting.c ${SOURCES})
target_link_libraries(bench_sorting PRIVATE OpenSSL::Crypto)
set_property(TARGET bench_sorting APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_1 -DBALANCED=1")
