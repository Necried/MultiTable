#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "rsqrt.h"
#include "float_defs.h"
#include "float_utils.h"
#include "utils.h"

// The Newton method for square root is
// y_n+1 = y_n * (3/2 - 1/2 * ay_n^2)
// The above equation computes sqrt
// To obtain sqrt: sqrt(a) = a * y
double calculate_ytest_dp_newton(unsigned __int128 a, __uint64_t *sqrtTable,
                       __uint64_t *linearTable, __uint8_t *aTable) {
  const unsigned __int128 leading_one = 1ull << fracBits;
  const __uint32_t FRACBITS_SP = 27;

  // The table lookup value is y0
  unsigned __int128 a_aligned = a >> (fracBits-FRACBITS_SP);
  unsigned __int128 x0 = three_table_procedure(a_aligned, sqrtTable, linearTable, aTable) << (fracBits-FRACBITS_SP);
  
  // iteration 1
  // Note that we have to be careful of the shifting here, as x0 * x0 uses
  // 2*fracBits bits
  unsigned __int128 x0_squared = big_mult(x0, x0);

  // Perform the rest of the computation with 2*fracBits,
  // which should give more than enough accuracy for the final result
  unsigned __int128 half_ax0_squared = (a * x0_squared) >> 1;
  unsigned __int128 three_over_two = (leading_one + (leading_one >> 1)) << fracBits;

  unsigned __int128 x1_rsqrt_inner = (three_over_two - half_ax0_squared) >> fracBits;
  unsigned __int128 x1_rsqrt = big_mult(x0, x1_rsqrt_inner);
  
  // iteration 2
  unsigned __int128 x1_squared = big_mult(x1_rsqrt, x1_rsqrt);

  unsigned __int128 half_ax1_squared = (a * x1_squared) >> 1;

  unsigned __int128 x2_rsqrt_inner = (three_over_two - half_ax1_squared) >> fracBits;
  unsigned __int128 x2_rsqrt = big_mult(x1_rsqrt, x2_rsqrt_inner);

  /*
  DEBUG_PRINT(("%-20s: %016llx %016llx\n",
	   "leading_one", FORMAT_UINT128(leading_one)));
  DEBUG_PRINT(("%-20s: %016llx %016llx %10.16f\n",
	   "a", FORMAT_UINT128(a), mantissa_to_double(a)));
  DEBUG_PRINT(("%-20s: %016llx %016llx %10.16f\n",
	   "x0", FORMAT_UINT128(x0), mantissa_to_double(x0)));
  DEBUG_PRINT(("%-20s: %016llx %016llx %10.16f\n",
	   "x0_squared", FORMAT_UINT128(x0_squared), mantissa_to_double(x0_squared)));
  DEBUG_PRINT(("%-20s: %016llx %016llx %10.16f\n",
	   "half_ax0_squared", FORMAT_UINT128(half_ax0_squared), mantissa_to_double(half_ax0_squared)));
  DEBUG_PRINT(("%-20s: %016llx %016llx %10.16f\n",
	   "three_over_two", FORMAT_UINT128(three_over_two), mantissa_to_double(three_over_two)));
  DEBUG_PRINT(("%-20s: %016llx %016llx %10.16f\n",
	   "x1_rsqrt_inner", FORMAT_UINT128(x1_rsqrt_inner), mantissa_to_double(x1_rsqrt_inner)));
  DEBUG_PRINT(("%-20s: %016llx %016llx %10.16f\n",
	   "x1_rsqrt", FORMAT_UINT128(x1_rsqrt), mantissa_to_double(x1_rsqrt)));
  DEBUG_PRINT(("%-20s: %016llx %016llx %10.16f\n",
	   "x2_rsqrt", FORMAT_UINT128(x1_rsqrt), mantissa_to_double(x2_rsqrt)));
  */
  return mantissa_to_double(x2_rsqrt);
}

double calculate_ytrue_dp(unsigned __int128 x){
  // printf("true: %lf\n", 1/(mantissa_to_double(x)));
  return 1/(sqrt(mantissa_to_double(x)));
}