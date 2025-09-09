#include "dtof_rounding.h"
#include <fenv.h>
#include <math.h>

float dtof_pinf(double x) {
  float res;
  int rndmode;
  // Save current rounding mode
  rndmode = fegetround();
  // Set rounding mode to round to nearest
  fesetround(FE_UPWARD);
  res = (float)x;
  // Restore rounding mode
  fesetround(rndmode);
  return res;
}
