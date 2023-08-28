#pragma once
#include <math.h>
#include <inttypes.h>

static inline
uint64_t x86_64_rtdsc(void) {
#if defined(__x86_64__) || defined(__x86_64) || defined(__amd64__) || defined(__amd64)
  unsigned long long result;
    __asm__ __volatile__(
        "rdtscp;"
        "shl $32, %%rdx;"
        "or %%rdx, %%rax"
        : "=a"(result)
        :
        : "%rcx", "%rdx");
  return result;
#else
#warning "rtdsc is not available on the target platform"
  return 0;
#endif
}

typedef struct {
     long double mean;
     long double M2;
     long count;
} welford_t;

static inline
void welford_init(welford_t* state){
    state->mean = 0.0;
    state->M2 = 0.0;
    state->count = 0;
    return;
}

static inline
void welford_update(welford_t* state, long double sample){
    long double delta, delta2;
    state->count = state->count + 1;
    delta = sample - state->mean;
    state->mean += delta / (long double)(state->count);
    delta2 = sample - state->mean;
    state->M2 += delta * delta2;
}

static inline
double welch_t_statistic(const welford_t state1,
                         const welford_t state2){
 long double num, den, var1, var2;
 var1 = state1.M2/(long double)(state1.count-1);
 var2 = state2.M2/(long double)(state2.count-1);

 num = state1.mean - state2.mean;
 den = sqrtl(var1/(long double) state1.count + var2/(long double) state2.count );

 return num/den;
}

static inline
void welford_print(const welford_t state){
     printf("%.2Lf,%.2Lf",
              state.mean,
              sqrtl(state.M2/(long double)(state.count-1)));
}

static inline
long double welford_stddev(const welford_t state){
     return sqrtl(state.M2/(long double)(state.count-1));
}

static inline
long double welford_mean(const welford_t state){
     return state.mean;
}
