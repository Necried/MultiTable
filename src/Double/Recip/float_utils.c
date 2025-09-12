#include "float_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "float_defs.h"

// Because the range is in [1,2), when mantissa is 0, assume the floating point double is 0;
double mantissa_to_double (unsigned __int128 mantissa) {

  //printf("mantissa %016llx %016llx \n",(unsigned long long)mantissa,(unsigned long long)condition);
  double result;
  result = ((double)mantissa)/pow(2.0,fracBits);

  //printf("mantissa %lf pow %lf\n", (double)mantissa, pow(2.0,fracBits));
  return result;
}

//mantissa with imply bit
//  double_to_mantissa works in a interval of [0.5,2.0)
unsigned __int128 double_to_mantissa (double x) {
  double_union_t ux;
  unsigned __int128 mantissa = 0;
  ux.d = x;
  // printf("exp %d %016llX ", fix.exp, fix.exp);
  int exp;
  exp = (ux.ull & 0x7FF0000000000000ull) >> 52;
  exp -= 1023; // subtract bias
  mantissa += (1ull << (fracBits+exp)) /*implybit*/ ;
  // int shift = x < 1.0 ?  FRACBITS - 1 : FRACBITS;
  mantissa += (fracBits+exp) > 52 ?
      (ux.ull & 0x000FFFFFFFFFFFFFull)<<((fracBits+exp)-52)
    : (ux.ull & 0x000FFFFFFFFFFFFFull)>>(52-(fracBits+exp));  // without assumed 1

  //printf("x %10.10f u %016llx exp %d %016x \n",x ,ux.ll,exp,exp);
  //fix = exp < 0 ? fix + (1ull << (FRACBITS-1)) : fix + (1ull << (FRACBITS)); // add imply 1
  //printf("fix %016llX \n",fix);
  return mantissa;
}

double mypow2 (int expo)
{
  double res;
  // Compute 2**expo being careful about pow() problems with denorms
  if (expo >= -1022) {
    res = pow (2.0, expo);
  }
  else {
    res = pow (2.0, expo+53) * pow (2.0, -53);
  }
  return res;
}

double ulp_d (double y)
{
  long long expobits;
  int expo;
  union {double d; unsigned long long ull;} u;
  u.d = y;
  // printf("u.ull    =%08llX\n", u.ull);

  expobits = u.ull & 0x7FF0000000000000ull;
  //printf("expobits=%08X\n", expobits);
  expo = expobits >> (64-12);
  //printf("expo    =%08X=%d\n", expo, expo);
  expo -= 1023;
  //printf("expo    =%08X=%d\n", expo, expo);
  if (expo >= -1022) {
    return mypow2 (expo-52);
  }
  else {
    printf("ulp_d denormalized\n");
    exit(1);
    return mypow2 (-149);
  }
}

float ulp_s (float y)
{
  unsigned int expobits;
  int expo;
  union {float f; int i;} u;
  u.f = y;
  //printf("u.i    =%08X\n", u.i);
  expobits = u.i & 0x7F800000;
  //printf("expobits=%08X\n", expobits);
  expo = expobits >> (32-9);
  //printf("expo    =%08X=%d\n", expo, expo);
  expo -= 127;
  //printf("expo    =%08X=%d\n", expo, expo);
  if (expo >= -126) {
    return mypow2 (expo-23);
  }
  else {
    return mypow2 (-149);
  }
}

double ieee_to_double(unsigned long long x) {
  union {unsigned long long l; double d;} tmp;
  tmp.l = x;
  return tmp.d;
}

unsigned __int128 big_mult(unsigned __int128 x, unsigned __int128 y) {
  const unsigned int halfFracBits = fracBits / 2;
  unsigned __int128 x0 = x & (1ull<<halfFracBits)-1;
  unsigned __int128 y0 = y & (1ull<<halfFracBits)-1;

  unsigned __int128 x1 = x >> halfFracBits;
  unsigned __int128 y1 = y >> halfFracBits;

  unsigned __int128 z0 = x0 * y0;
  z0 >>= fracBits;
  unsigned __int128 z1 = x1 * y0 + x0 * y1;
  z1 >>= halfFracBits;
  unsigned __int128 z2 = x1 * y1;

  return z2 + z1 + z0;
}
