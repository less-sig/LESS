#
# LESS Profiling makefile
#

# user definable variables

# NONE

##########################################
# C Compiliation
CC ?= /usr/bin/cc
# CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
# 	-Wshadow -Wpointer-arith -O3 -fomit-frame-pointer 

CFLAGS += -O3 
# LD_FLAGS+= -Wl,--stack-size=1000000

##########################################
# C Sources
LESS_SRC_DIR    := lib
LESS_HEADER_DIR := include

INCLUDE_DIR  = -I$(LESS_HEADER_DIR)

SOURCES_BENCH = lib/test/less_test.c

SOURCES_LESS = $(LESS_SRC_DIR)/codes.c \
				$(LESS_SRC_DIR)/utils.c \
				$(LESS_SRC_DIR)/fips202.c \
				$(LESS_SRC_DIR)/keccakf1600.c \
				$(LESS_SRC_DIR)/LESS.c \
				$(LESS_SRC_DIR)/monomial.c \
				$(LESS_SRC_DIR)/rng.c \
				$(LESS_SRC_DIR)/seedtree.c \
				$(LESS_SRC_DIR)/sign.c 

HEADERS_LESS = $(LESS_HEADER_DIR)/api.h \
				$(LESS_HEADER_DIR)/parameters.h \
				$(LESS_HEADER_DIR)/utils.h \
				$(LESS_HEADER_DIR)/codes.h \
				$(LESS_HEADER_DIR)/fips202.h \
				$(LESS_HEADER_DIR)/fq_arith.h \
				$(LESS_HEADER_DIR)/keccakf1600.h \
				$(LESS_HEADER_DIR)/LESS.h \
				$(LESS_HEADER_DIR)/monomial_mat.h \
				$(LESS_HEADER_DIR)/rng.h \
				$(LESS_HEADER_DIR)/seedtree.h \
				$(LESS_HEADER_DIR)/sha3.h 

all: less1b less1i less1s less3b less3s less5b less5s

bench: $(SOURCES_LESS) $(HEADERS_LESS) $(SOURCES_BENCH)
	$(CC) -pg $(CFLAGS) $(LD_FLAGS) -DCATEGORY_1=1 -DBALANCED=1  $(INCLUDE_DIR) $(SOURCES_LESS) $(SOURCES_BENCH) -o test

# clean:
# 	rm -rf *.so correction