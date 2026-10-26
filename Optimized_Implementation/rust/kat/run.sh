#!/bin/sh
# Builds the C reference implementation together with `less_kat.c` and prints
# the reference vectors on stdout.
#
# NOTE: nothing outside this directory is written to, the object files and the
# binary go into a temporary build directory.
#
# usage: kat/run.sh [CATEGORY] [TARGET]      defaults to 252 192
set -e

CATEGORY=${1:-252}
TARGET=${2:-192}

HERE=$(cd "$(dirname "$0")" && pwd)
REF=$HERE/../../../Reference_Implementation
OUT=${KAT_BUILD_DIR:-$(mktemp -d)}

cc -O2 -std=c99 -w \
   -DCATEGORY="$CATEGORY" -DTARGET="$TARGET" \
   -I"$REF/include" \
   "$HERE/less_kat.c" \
   "$REF/lib/codes.c" \
   "$REF/lib/fips202.c" \
   "$REF/lib/keccakf1600.c" \
   "$REF/lib/LESS.c" \
   "$REF/lib/monomial.c" \
   "$REF/lib/rng.c" \
   "$REF/lib/seedtree.c" \
   "$REF/lib/sign.c" \
   "$REF/lib/utils.c" \
   "$REF/lib/sort.c" \
   "$REF/lib/canonical.c" \
   "$REF/lib/transpose.c" \
   -o "$OUT/less_kat"

"$OUT/less_kat"
