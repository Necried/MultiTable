#include "dtof_rounding.h"
#include <fenv.h>
#include <math.h>

float dtof_zero(double x) {
  float res;
  int rndmode;
  // Save current rounding mode
  rndmode = fegetround();
  // Set rounding mode to round to nearest
  fesetround(FE_TOWARDZERO);
  res = (float)x;
  // Restore rounding mode
  fesetround(rndmode);
  return res;
}
