cmake_minimum_required(VERSION 3.22)
project(less_avx2)
enable_language(C ASM)
set(CMAKE_C_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE ON)

# needed for nvim
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
if(CMAKE_EXPORT_COMPILE_COMMANDS)
  set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES
          ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()

set(ALLOWED_WARNINGS      " -Wno-unused-function -Wno-ignored-attributes -Wno-type-limits -Wno-gnu-binary-literal -fno-inline")
set(CMAKE_C_FLAGS_DEBUG   " -DUSE_AVX2 -g -O0 -DDEBUG  -Wall -Wextra -Wpedantic -mavx2 -mavx -march=native ${ALLOWED_WARNINGS}")
set(CMAKE_C_FLAGS_RELEASE " -DUSE_AVX2 -g -O3 -DNDEBUG -Wall -Wextra -Wpedantic -mavx2 -mavx -march=native -flto -ftree-vectorize -funroll-loops -fomit-frame-pointer ${ALLOWED_WARNINGS}")


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
        ${PROJECT_SOURCE_DIR}/lib/canonical.c
        ${PROJECT_SOURCE_DIR}/lib/sort.c
        ${PROJECT_SOURCE_DIR}/lib/transpose.c
        ${PROJECT_SOURCE_DIR}/lib/KeccakP-1600-AVX2.s
        ${PROJECT_SOURCE_DIR}/lib/transpose_avx2.c
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
        ${PROJECT_SOURCE_DIR}/include/canonical.h
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
  foreach(optimiz_target ${PARAM_TARGETS})
    # settings for benchmarking binary
    set(TARGET_BINARY_NAME LESS_benchmark_cat_${category}_${optimiz_target})
    add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES} ${PROJECT_SOURCE_DIR}/lib/bench/less_benchmark.c)
    target_link_libraries(${TARGET_BINARY_NAME} m)
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_${category}=1 -D${optimiz_target}=1 ${KECCAK_EXTERNAL_ENABLE} ")

    # settings for unit tests binary
    set(TARGET_BINARY_NAME LESS_test_cat_${category}_${optimiz_target})
    add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES} ${PROJECT_SOURCE_DIR}/lib/test/less_test.c)
    target_link_libraries(${TARGET_BINARY_NAME} m)
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_${category} -D${optimiz_target}=1 ${KECCAK_EXTERNAL_ENABLE} ")

    # KATS generation
    set(TARGET_BINARY_NAME LESS_nist_cat_${category}_${optimiz_target})
    add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES}  ${PROJECT_SOURCE_DIR}/lib/nist/KAT_NIST_rng.c ${PROJECT_SOURCE_DIR}/lib/nist/PQCgenKAT_sign.c)
    target_link_libraries(${TARGET_BINARY_NAME} m crypto)
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_${category} -D${optimiz_target}=1 ${KECCAK_EXTERNAL_ENABLE} ")

    # TODO remove once everything is done
    set(TARGET_BINARY_NAME LESS_test_cf_${category}_${optimiz_target})
    add_executable(${TARGET_BINARY_NAME} ${PROJECT_SOURCE_DIR}/lib/test/test_canonical.c ${SOURCES})
    target_link_libraries(${TARGET_BINARY_NAME} m crypto)
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_${category} -D${optimiz_target}=1 ${KECCAK_EXTERNAL_ENABLE} ")

    # TODO remove once everything is done
    set(TARGET_BINARY_NAME LESS_bench_cf_${category}_${optimiz_target})
    add_executable(${TARGET_BINARY_NAME} ${PROJECT_SOURCE_DIR}/lib/bench/bench_cf.c ${SOURCES})
    target_link_libraries(${TARGET_BINARY_NAME} m crypto)
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_${category} -D${optimiz_target}=1 ${KECCAK_EXTERNAL_ENABLE} ")

  endforeach(optimiz_target)

  # TODO remove once everything is done
  set(TARGET_BINARY_NAME LESS_test_arith_${category})
  add_executable(${TARGET_BINARY_NAME} ${PROJECT_SOURCE_DIR}/lib/test/test_arith.c ${SOURCES})
  target_link_libraries(${TARGET_BINARY_NAME} m crypto)
  set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_${category} -DSHORT_SIG=1 ${KECCAK_EXTERNAL_ENABLE} ")

  # TODO remove once everything is done
  set(TARGET_BINARY_NAME LESS_test_sort_${category})
  add_executable(${TARGET_BINARY_NAME} ${PROJECT_SOURCE_DIR}/lib/test/sort.c ${SOURCES})
  target_link_libraries(${TARGET_BINARY_NAME} m crypto)
  set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_${category} -DSHORT_SIG=1 ${KECCAK_EXTERNAL_ENABLE} ")

  # TODO remove once everything is done
  set(TARGET_BINARY_NAME LESS_bench_sort_${category})
  add_executable(${TARGET_BINARY_NAME} ${PROJECT_SOURCE_DIR}/lib/bench/bench_sorting.c ${SOURCES})
  target_link_libraries(${TARGET_BINARY_NAME} m crypto)
  set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY_${category} -DSHORT_SIG=1 ${KECCAK_EXTERNAL_ENABLE} ")
endforeach(category)

