#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "recip.h"
#include "float_defs.h"
#include "float_utils.h"
#include "utils.h"

// Source: Series Approximation Methods for Divide and Square Root in the Power3 Processor
// Authors: Agarwal, Gustavson, Schmookler, IBM Corporation
// Idea: A modified Taylor's series third-order polynomial, with
// intermediate computations only requiring 1 FMA and 1 FM for the polynomial computation
// However, we can accelerate the computations with the fact that the
// initial residual (e) and residual squared (e^2) are very small, compared
// to the original paper, due to the nature of our tables being more accurate
// Refer to p.5, Table 1, SP Divide in the paper for the original implementation.
// We conjecture that DP Divide is too much accuracy for what we have.
// Inputs: Denominator b, tables
double calculate_ytest_dp(unsigned __int128 b, __uint64_t *recipTable,
                       __uint64_t *linearTable, __uint8_t *aTable) {
  const unsigned __int128 leading_one = 1ull << fracBits;
  const __uint32_t FRACBITS_SP = 31;

  // The table lookup value is y0
  unsigned __int128 b_aligned = b >> (fracBits-FRACBITS_SP);
  unsigned __int128 y0 = three_table_procedure(b_aligned, recipTable, linearTable, aTable) << (fracBits-FRACBITS_SP);

  // First we compute the residual e, and residual_squared e^2
  // Extra step: We do 1 - b*y0 if b*y0 < 1, else do b*y0 - 1
  unsigned __int128 r = big_mult(b, y0);
  bool r_lesser_than_e = r < leading_one;
  unsigned __int128 e = 
       r_lesser_than_e ? leading_one - r : r - leading_one;
  unsigned __int128 e_squared = big_mult(e, e);

  // Calculate t1 and t3.
  unsigned __int128 t1 = leading_one + (r_lesser_than_e ? e : -e) + e_squared; // leading_one + e_squared;
  unsigned __int128 t3 = y0; // r_lesser_than_e ? y0 + big_mult(y0, e) : y0 - big_mult(y0, e);

  // The final result q1
  unsigned __int128 q1 = big_mult(t1, t3);

  // Print debugs
  DEBUG_PRINT(("%-20s: %016llx %016llx\n",
	   "leading_one", FORMAT_UINT128(leading_one)));
  DEBUG_PRINT(("%-20s: %016llx %016llx %10.16f\n%-20s: %016llx %016llx\n",
	   "y0", FORMAT_UINT128(y0), mantissa_to_double(y0),
	   "y0 (shifted)", FORMAT_UINT128(y0<<fracBits)));
  DEBUG_PRINT(("%-20s: %016llx %016llx\n%-20s: %016llx %016llx\n%-20s: %016llx %016llx\n",
	   "b*y0", FORMAT_UINT128(r),
	   "e = 1 - b*y0", FORMAT_UINT128(e),
	   "e_squared", FORMAT_UINT128(e_squared)));
  DEBUG_PRINT(("%-20s: %016llx %016llx\n%-20s: %016llx %016llx\n",
	   "t1 = 1 + e^2", FORMAT_UINT128(t1),
	   "t3 = y0 + y0*e", FORMAT_UINT128(t3)));
  DEBUG_PRINT(("%-20s: %016llx %016llx\n",
	   "q1 = t1 * t3", FORMAT_UINT128(q1)));

  // Normalize and return the result
  return mantissa_to_double(q1);
}

double calculate_ytest_dp_newton(unsigned __int128 a, __uint64_t *recipTable,
                       __uint64_t *linearTable, __uint8_t *aTable) {
  const unsigned __int128 leading_one = 1ull << fracBits;

  // The table lookup value is y0
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