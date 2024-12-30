#ifndef CYCLES_H
#define CYCLES_H

#include <stdint.h>

#if defined(__APPLE__) && (defined(__aarch64__) || defined(__arm64__))
#define MACOS_KPERF
#include "m1cycles.h"
#endif


static
void setup_cycle_counter(void)
{
#ifdef MACOS_KPERF
   kperf_init_once();
#endif
}

static inline
uint64_t read_cycle_counter(void)
{
#ifdef MACOS_KPERF
   return __m1_rdtsc();
#elif defined(__amd64) || defined(__amd64__) || defined(__x86_64__)

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
  return 0;
#endif
}

#endif //CYCLES_H
