cmake_minimum_required(VERSION 3.5)
project(mirith)
enable_language(C ASM)
set(CMAKE_C_STANDARD 11)

set(CMAKE_VERBOSE_MAKEFILE ON)

# needed for nvim
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
if(CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES
            ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()

if(APPLE)
set(CMAKE_C_FLAGS_DEBUG "-g -DDEBUG -O0 -Wall -flax-vector-conversions -fsanitize=address")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -Wall -flto -march=native -flax-vector-conversions -ftree-vectorize -funroll-loops -fomit-frame-pointer -fno-stack-protector")
else()
set(CMAKE_C_FLAGS_DEBUG "-g -DDEBUG -O0 -Wall -fsanitize=address")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -Wall -flto -march=native -ftree-vectorize -funroll-loops -fomit-frame-pointer -fno-stack-protector")
endif()

include_directories(lib/less)


set(SOURCES
        ${PROJECT_SOURCE_DIR}/lib/less/codes.c
        ${PROJECT_SOURCE_DIR}/lib/less/fips202.c
        ${PROJECT_SOURCE_DIR}/lib/less/keccakf1600.c
        ${PROJECT_SOURCE_DIR}/lib/less/LESS.c
        ${PROJECT_SOURCE_DIR}/lib/less/monomial.c
        ${PROJECT_SOURCE_DIR}/lib/less/rng.c
        ${PROJECT_SOURCE_DIR}/lib/less/seedtree.c
        ${PROJECT_SOURCE_DIR}/lib/less/sign.c
        ${PROJECT_SOURCE_DIR}/lib/less/utils.c
        ${PROJECT_SOURCE_DIR}/lib/less/sort.c
        ${PROJECT_SOURCE_DIR}/lib/less/canonical.c
        ${PROJECT_SOURCE_DIR}/lib/less/transpose.c
)

set(HEADERS
        ${PROJECT_SOURCE_DIR}/lib/less/api.h
        ${PROJECT_SOURCE_DIR}/lib/less/codes.h
        ${PROJECT_SOURCE_DIR}/lib/less/fips202.h
        ${PROJECT_SOURCE_DIR}/lib/less/keccakf1600.h
        ${PROJECT_SOURCE_DIR}/lib/less/sha3.h
        ${PROJECT_SOURCE_DIR}/lib/less/fq_arith.h
        ${PROJECT_SOURCE_DIR}/lib/less/LESS.h
        ${PROJECT_SOURCE_DIR}/lib/less/monomial_mat.h
        ${PROJECT_SOURCE_DIR}/lib/less/parameters.h
        ${PROJECT_SOURCE_DIR}/lib/less/rng.h
        ${PROJECT_SOURCE_DIR}/lib/less/seedtree.h
        ${PROJECT_SOURCE_DIR}/lib/less/utils.h
        ${PROJECT_SOURCE_DIR}/lib/less/sort.h
        ${PROJECT_SOURCE_DIR}/lib/less/transpose.h
)

include_directories(lib/less)

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
    set(TARGET_BINARY_NAME LESS_benchmark_${category}_${optimize_target})
    add_executable(${TARGET_BINARY_NAME} ${HEADERS} ${SOURCES} ${PROJECT_SOURCE_DIR}/src/main.c)
    target_link_libraries(${TARGET_BINARY_NAME} m)
    set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY COMPILE_FLAGS "-DCATEGORY=${category} -DTARGET=${optimize_target} ${KECCAK_EXTERNAL_ENABLE} ")
    target_include_directories(${TARGET_BINARY_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lib/test)
  endforeach(optimize_target)
endforeach(category)
