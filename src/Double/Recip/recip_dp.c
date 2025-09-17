#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "recip.h"
#include "float_defs.h"
#include "float_utils.h"
#include "utils.h"

// The Newton method for recip is
// y_n+1 = y_n * (2 - a*y_n)
// Calculates two iteration of Newton-Rhapson after the initial three-table method.
double calculate_ytest_dp_newton(unsigned __int128 a, __uint64_t *recipTable,
                       __uint64_t *linearTable, __uint8_t *aTable) {
  // The table lookup value is x0
  unsigned __int128 a_aligned = a >> (fracBits-31);
  unsigned __int128 x0 = three_table_procedure(a_aligned, recipTable, linearTable, aTable) << (fracBits-31);

  // iteration 1
  // Note that we have to be careful of the shifting here, as a * x0 uses
  // 2*fracBits bits
  unsigned __int128 ax0 = a * x0;

  // WARNING: We cannot do 2ull<<2*fracBits, as unsigned __128 doesn't
  // shift correctly past 64 bits... so it has to be split
  unsigned __int128 two_sub_ax0 = 2ull<<fracBits; // (2 - a*x0)
  two_sub_ax0 <<= fracBits;
  two_sub_ax0 -= ax0;
  two_sub_ax0 >>= fracBits;

  unsigned __int128 x1 = big_mult(x0, two_sub_ax0); // x1 = x0 * (2 - a*x0)

  // iteration 2
  unsigned __int128 ax1 = a * x1;

  unsigned __int128 two_sub_ax1 = 2ull<<fracBits; // (2 - a*x1)
  two_sub_ax1 <<= fracBits;
  two_sub_ax1 -= ax1;
  two_sub_ax1 >>= fracBits;

  unsigned __int128 x2 = big_mult(x1, two_sub_ax1); // x2 = x1 * (2 - a*x1)

  return mantissa_to_double(x2);
}

double calculate_ytrue_dp(unsigned __int128 x){
  // printf("true: %lf\n", 1/(mantissa_to_double(x)));
  return 1/(mantissa_to_double(x));
}