/* TODO: convert ints to longs for base counts and such things? */

#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include <assert.h>
#include <math.h>

#include "complexity.h"

#define RND MPFR_RNDN
#define NUM_BASES 4

/* calculates max! and stores the result in result */
void fact(mpfr_t result, unsigned long long max)
{
  /* need to check that it is initialized */

  /* set the result to 1 to start */
  mpfr_set_ui(result, 1, RND);

  unsigned long long i = 0;
  for (i = 2; i <= max; i++) {
    /* result = result * i */
    mpfr_mul_ui(result, result, i, RND);
  }
}

/* log base 4 of x -> log(x) / log(4), sets result to the result */
void log4(mpfr_t result, mpfr_t x) {
  /* make sure values are initialized */

  mpfr_t val1;
  mpfr_init(val1);
  mpfr_set_d(val1, 1, RND);

  mpfr_t val2;
  mpfr_init(val2);
  mpfr_set_d(val2, 1, RND);

  mpfr_t base;
  mpfr_init(base);
  mpfr_set_ui(base, 4, RND);

  mpfr_set_ui(result, 1, RND);


  /* set val1 and val2 to the proper logs */
  mpfr_log(val1, x, RND);
  mpfr_log(val2, base, RND);

  assert(mpfr_cmp_d(val2, 0.0) != 0); // can the log ever be 0?

  mpfr_div(result, val1, val2, RND);

  mpfr_clear(val1);
  mpfr_clear(val2);
  mpfr_clear(base);
}

double cwf(unsigned long long base_counts[], unsigned long long window)
{
  assert(window != 0);
  assert(base_counts != NULL);

  unsigned long long i = 0;
  unsigned long long count = 0;

  mpfr_t val;
  mpfr_init(val);
  mpfr_set_ui(val, 1, RND);

  mpfr_t loop_result;
  mpfr_init(loop_result);
  mpfr_set_ui(loop_result, 1, RND);

  mpfr_t log_result;
  mpfr_init(log_result);
  mpfr_set_ui(log_result, 1, RND);


  for (i = 0; i < NUM_BASES; i++) {
    count = base_counts[i];
    assert(count > 0);

    fact(val, count);

    /* loop_result *= val */
    mpfr_mul(loop_result, loop_result, val, RND);
  }

  assert(mpfr_cmp_ui(loop_result, 0) != 0);

  /* val = factorial(window) */
  fact(val, window);
  /* val = val / loop_result */
  mpfr_div(val, val, loop_result, RND);

  /* log_result = log4(val) */
  log4(log_result, val);

  /* val = log_result * (1.0 / window) */
  mpfr_mul_d(val, log_result, (1.0 / window), RND);

  double output = 0.0;
  output = mpfr_get_d(val, RND);

  mpfr_clear(val);
  mpfr_clear(loop_result);
  mpfr_clear(log_result);

  return output;
}

double ce(unsigned long long base_counts[], unsigned long long window)
{
  assert(window != 0);
  unsigned long long i = 0;
  unsigned long long count = 0;
  double val = 0;
  double accum = 0;

  for (i = 0; i < NUM_BASES; i++) {
    count = base_counts[i];
    assert(count > 0);

    val = count / (double) window;
    accum += val * log(val) / log(4);
  }

  return -accum;
}

int main()
{
  unsigned long long base_counts[] = { 5, 5, 5, 5 };
  unsigned long long window = 20;
  double val = 0.0;

  val = cwf(base_counts, window);
  printf("%f\n", val);

  val = ce(base_counts, window);
  printf("%f\n", val);

  mpfr_free_cache();

  return 0;
}
