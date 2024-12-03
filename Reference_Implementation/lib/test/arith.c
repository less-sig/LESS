#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "../bench/cycles.h"

#define Q 127
#define Qm1 (126)
#define M 0x204081020408103ull

#define FQ_ELEM uint8_t
#define FQ_DOUBLEPREC uint16_t

static inline
FQ_ELEM fq_add1(const FQ_ELEM x, const FQ_ELEM y) {
   return (x + y) % Q;
}
static inline FQ_ELEM fq_mul1(const FQ_ELEM x,
                             const FQ_ELEM y) {
   return ((FQ_DOUBLEPREC)x * (FQ_DOUBLEPREC)y) % Q;
}

static inline FQ_ELEM fq_red1(FQ_DOUBLEPREC x) {
   return ((FQ_DOUBLEPREC) Q+x) % (FQ_DOUBLEPREC) Q;
}

static inline
FQ_ELEM fq_mul(FQ_ELEM a, FQ_ELEM b){
   FQ_DOUBLEPREC lo, hi;
   lo = a * b;
   hi = lo >> 7;
   lo += hi;
   hi <<= 7;
   return (lo - hi);
}

/*
 * Barrett reduction for uint8_t with prime Q = 127
 */
static inline
FQ_ELEM fq_red(FQ_ELEM a) {
   FQ_ELEM t;
   t = a >> 7;
   t &= Q;
   a += t;
   a &= Q;
   return a;
}

static inline
FQ_ELEM fq_add(const FQ_ELEM x, const FQ_ELEM y) {
   return fq_red(x + y);
}

/// super generic a \in [0, 2**16)
uint8_t fq_red2(const uint16_t a) {
   uint64_t lowbits = M * a;
   uint32_t highbits = ((__uint128_t)lowbits * Q) >> 64;
   return highbits - (Qm1 & (a >> 31));
}

#define ITERS (1u << 20)
int bench_sorting(void) {
   uint8_t a = 1, b = 2, c = 3, d;
    uint64_t cyc = 0, c1;
    for (uint64_t i = 0; i < ITERS; i++) {

        cyc -= x86_64_rtdsc();
        a = fq_red2(c + b);
        b = fq_red2(c + a);
        c = fq_red2(a + b);
        a = fq_red2(c + b);
        b = fq_red2(c + a);
        c = fq_red2(a + b);
        a = fq_red2(c + b);
        b = fq_red2(c + a);
        c = fq_red2(a + b);
        a = fq_red2(c + b);
        b = fq_red2(c + a);
        c = fq_red2(a + b);
        a = fq_red2(c + b);
        b = fq_red2(c + a);
        c = fq_red2(a + b);
        a = fq_red2(c + b);
        b = fq_red2(c + a);
        c = fq_red2(a + b);
        a = fq_red2(c + b);
        b = fq_red2(c + a);
        c = fq_red2(a + b);
        a = fq_red2(c + b);
        b = fq_red2(c + a);
        c = fq_red2(a + b);
        cyc += x86_64_rtdsc();
    }

    d = c;
    c1 = cyc/ITERS;
    printf("red2: %ld cyc\n", c1);

    cyc = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        cyc -= x86_64_rtdsc();
        a = fq_red1(c + b);
        b = fq_red1(c + a);
        c = fq_red1(a + b);
        a = fq_red1(c + b);
        b = fq_red1(c + a);
        c = fq_red1(a + b);
        a = fq_red1(c + b);
        b = fq_red1(c + a);
        c = fq_red1(a + b);
        a = fq_red1(c + b);
        b = fq_red1(c + a);
        c = fq_red1(a + b);
        a = fq_red1(c + b);
        b = fq_red1(c + a);
        c = fq_red1(a + b);
        a = fq_red1(c + b);
        b = fq_red1(c + a);
        c = fq_red1(a + b);
        a = fq_red1(c + b);
        b = fq_red1(c + a);
        c = fq_red1(a + b);
        a = fq_red1(c + b);
        b = fq_red1(c + a);
        c = fq_red1(a + b);
        cyc += x86_64_rtdsc();
    }

    d += c;
    cyc = cyc/ITERS;
    printf("red1: %ld cyc\n", cyc);
    printf("factor %lf\n", (double)cyc/(double)c1);
   return d;
}

int main(void) {
   // uint64_t m = computeM_s32(127);
   return bench_sorting();

   for (uint32_t i = 0; i < (127*127); i++) {
         const FQ_ELEM c2 = fq_red2(i);
         const FQ_ELEM c1 = fq_red1(i);
         assert(c1 == c2);
   }

   for (uint32_t i = 1; i < 128; i++) {
      for (uint32_t j = 127; j < 128; j++) {
         const FQ_ELEM c1 = fq_mul1(i, j);
         const FQ_ELEM c2 = fq_mul(i, j);
         assert(c1 == c2);
      }
   }
}
