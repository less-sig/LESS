cmake_minimum_required(VERSION 3.22)
project(less_neon)
enable_language(C ASM)
set(CMAKE_C_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE ON)

if(APPLE)
    include_directories("/opt/homebrew/include/")
    link_directories("/opt/homebrew/lib")
endif()

enable_testing()
find_package(OpenSSL)

# needed for nvim
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
if(CMAKE_EXPORT_COMPILE_COMMANDS)
  set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES
          ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()

set(ALLOWED_WARNINGS      " -Wno-unused-function -Wno-ignored-attributes -Wno-type-limits -Wno-gnu-binary-literal")
set(CMAKE_C_FLAGS_DEBUG "-g -DUSE_NEON -DDEBUG -O0  -march=native -flax-vector-conversions")
set(CMAKE_C_FLAGS_RELEASE "-g -DUSE_NEON -O3 -DNDEBUG -Wall -Wextra -Wpedantic -fpermissive -march=native -flax-vector-conversions -ftree-vectorize -funroll-loops -fomit-frame-pointer -fno-stack-protector ${ALLOWED_WARNING}")

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
        ${PROJECT_SOURCE_DIR}/lib/transpose_neon.c
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
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY=${category} -DTARGET=${optimize_target} ${KECCAK_EXTERNAL_ENABLE} ")
    target_include_directories(${TARGET_BINARY_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lib/test)

    # settings for unit tests binary
    set(TARGET_BINARY_NAME LESS_test_cat_${category}_${optimize_target})
    add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES} ${PROJECT_SOURCE_DIR}/lib/test/less_test.c)
    target_link_libraries(${TARGET_BINARY_NAME} m)
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY=${category} -DTARGET=${optimize_target} ${KECCAK_EXTERNAL_ENABLE} ")
    add_test(${TARGET_BINARY_NAME} ${TARGET_BINARY_NAME})
    target_include_directories(${TARGET_BINARY_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lib/test)

    # KATS generation
    set(TARGET_BINARY_NAME LESS_nist_cat_${category}_${optimize_target})
    add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES}  ${PROJECT_SOURCE_DIR}/lib/nist/KAT_NIST_rng.c ${PROJECT_SOURCE_DIR}/lib/nist/PQCgenKAT_sign.c)
    target_link_libraries(${TARGET_BINARY_NAME} PRIVATE crypto)
    include_directories(${OPENSSL_INCLUDE_DIR})
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY=${category} -DTARGET=${optimize_target} ${KECCAK_EXTERNAL_ENABLE} ")
    add_test(${TARGET_BINARY_NAME} ${TARGET_BINARY_NAME})
  endforeach(optimize_target)

  set(TARGET_BINARY_NAME LESS_test_arith_${category})
  add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES}  ${PROJECT_SOURCE_DIR}/lib/test/test_arith.c)
  include_directories(${OPENSSL_INCLUDE_DIR})
  set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY=${category} -DTARGET=192 ${KECCAK_EXTERNAL_ENABLE} ")
  add_test(${TARGET_BINARY_NAME} ${TARGET_BINARY_NAME})
endforeach(category)

