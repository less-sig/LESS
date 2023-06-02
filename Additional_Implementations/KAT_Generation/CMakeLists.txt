cmake_minimum_required(VERSION 3.7)

project(LESS C)

set(CC gcc)
# set(CC clang)

set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -Wall -pedantic -Wuninitialized -march=native -O3 -g3")
set(SANITIZE "")

message("Compilation flags:" ${CMAKE_C_FLAGS})

set(CMAKE_C_STANDARD 11)

find_library(KECCAK_LIB keccak)
if(NOT KECCAK_LIB)
 set(STANDALONE_KECCAK 1)
endif()

# default compilation picks reference codebase
set(REFERENCE 1)

# selection of specialized compilation units differing between ref and opt
# implementations. Sample codes.c file being used as a placekeeper for now
if(REFERENCE)
set(BASE_DIR ../../Reference_Implementation)
set(SPEC_HEADERS
        ${BASE_DIR}/lib/codes.c
)
else()
set(BASE_DIR ../../Optimized_Implementation)
set(SPEC_HEADERS
        ${BASE_DIR}/lib/codes.c
)
endif()

set(HEADERS
    ${SPEC_HEADERS}
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
    ${BASE_DIR}/lib/LESS.c
    ${BASE_DIR}/lib/monomial.c
    ${BASE_DIR}/lib/rng.c
    ${BASE_DIR}/lib/seedtree.c
    ${BASE_DIR}/lib/utils.c
    ${BASE_DIR}/lib/sign.c
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


foreach(category RANGE 1 5 2)
    if(category EQUAL 1)
    set(PARAM_TARGETS SIG_SIZE BALANCED PK_SIZE)
    else()
    set(PARAM_TARGETS SIG_SIZE PK_SIZE)
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
                                    ./include)
         target_link_libraries(${TARGET_BINARY_NAME} m ssl crypto ${SANITIZE} ${KECCAK_EXTERNAL_LIB})
         set_target_properties(${TARGET_BINARY_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./bin)
         set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY
            COMPILE_FLAGS "-DCATEGORY_${category} -D${optimiz_target}=1 ${KECCAK_EXTERNAL_ENABLE} ")
    endforeach(optimiz_target)
endforeach(category)
