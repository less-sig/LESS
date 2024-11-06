#ifndef CYCLES_H
#define CYCLES_H

#if defined(__aarch64__) || defined(_M_ARM64)
#ifdef __APPLE__
#define _MAC_OS_
#define _M1CYCLES_
#include "m1cycles.h"
#endif

// TODO rename
static inline
uint64_t x86_64_rtdsc(void)
{
   unsigned long long result = __m1_rdtsc();;
   return result;
}
#else
static inline
uint64_t x86_64_rtdsc(void) {
    unsigned long long result;
    __asm__ __volatile__(
            "rdtscp;"
            "shl $32, %%rdx;"
            "or %%rdx, %%rax"
            : "=a"(result)
            :
            : "%rcx", "%rdx");
    return result;
}
#endif

#endif //CYCLES_H
