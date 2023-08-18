/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
* @version 1.1 (March 2023)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
 *
 * This code is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 **/

#pragma once
#include <math.h>
#include <stdint.h>
static inline
uint64_t x86_64_rtdsc(void)
{
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

typedef struct {
   long double mean;
   long double M2;
   long count;
} welford_t;

static inline
void welford_init(welford_t *state)
{
   state->mean = 0.0;
   state->M2 = 0.0;
   state->count = 0;
   return;
}

static inline
void welford_update(welford_t *state, long double sample)
{
   long double delta, delta2;
   state->count = state->count + 1;
   delta = sample - state->mean;
   state->mean += delta / (long double)(state->count);
   delta2 = sample - state->mean;
   state->M2 += delta * delta2;
}

static inline
double welch_t_statistic(const welford_t state1,
                         const welford_t state2)
{
   long double num, den, var1, var2;
   var1 = state1.M2/(long double)(state1.count-1);
   var2 = state2.M2/(long double)(state2.count-1);

   num = state1.mean - state2.mean;
   den = sqrtl(var1/(long double) state1.count + var2/(long double) state2.count );

   return num/den;
}

static inline
void welford_print(const welford_t state)
{
   printf("%.2Lf,%.2Lf",
          state.mean,
          sqrtl(state.M2/(long double)(state.count-1)));
}

static inline
long double welford_stddev(const welford_t state)
{
   return sqrtl(state.M2/(long double)(state.count-1));
}

static inline
long double welford_mean(const welford_t state)
{
   return state.mean;
}
