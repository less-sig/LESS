#ifndef CYCLES_H
#define CYCLES_H

#include <stdint.h>

#if defined(__APPLE__) && (defined(__aarch64__) || defined(__arm64__))
#include "m1cycles.h"
#endif

__attribute__((unused))
static
void setup_cycle_counter(void)
{
#if defined(__APPLE__) && (defined(__aarch64__) || defined(__arm64__))
    __m1_setup_rdtsc();
#endif
}

__attribute__((unused))
static inline
uint64_t read_cycle_counter(void)
{
#if defined(__APPLE__) && (defined(__aarch64__) || defined(__arm64__))
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
