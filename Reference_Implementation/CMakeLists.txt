cmake_minimum_required(VERSION 3.22)
project(less)
enable_language(C ASM)
set(CMAKE_C_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE ON)

# needed for nvim
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
if(CMAKE_EXPORT_COMPILE_COMMANDS)
  set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES
          ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()

set(ALLOWED_WARNINGS " -Wno-unused-function -Wno-type-limits")

IF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm") 
	set(CMAKE_C_FLAGS_DEBUG "-g -DDEBUG -O0  -march=native -flax-vector-conversions ${ALLOWED_WARNINGS}")
	set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -Wall -Wextra -Wpedantic -fpermissive -march=native -flax-vector-conversions -ftree-vectorize -funroll-loops -fomit-frame-pointer -fno-stack-protector ${ALLOWED_WARNINGS}")
else()
	set(CMAKE_C_FLAGS_DEBUG "-g -DDEBUG -O0 -Wall -fpermissive  -Wpedantic -Wextra -fpermissive -mavx2 -mavx -fsanitize=address ${ALLOWED_WARNINGS}")
	set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -Wall -Wextra -Wpedantic -fpermissive -march=native -mavx2 -mavx -ftree-vectorize -funroll-loops -fomit-frame-pointer -fno-stack-protector ${ALLOWED_WARNINGS}")
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
  endforeach(optimiz_target)
endforeach(category)
