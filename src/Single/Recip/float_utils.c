#include "float_utils.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

double float_internal_to_double(float_internal fix) {
  double x;
  if (fix.exp >= -1022) {                          // normalized
    x = fix.significand * pow(2.0, -1 * FRACBITS); // fraction
    x += 1;                                        // assumed 1
    x *= pow(2.0, fix.exp);                        // exponent
  } else {                                         // denormalized
    printf("float_internal_to_double:  denormalized x not implemented\n");
    assert(0);
  }
  return x;
}

// Because the range is in [1,2), when mantissa is 0, assume the floating point
// double is 0;
double mantissa_to_double(unsigned __int128 mantissa) {

  // printf("mantissa %016llx %016llx \n",(unsigned long long)mantissa,(unsigned
  // long long)condition);
  double result;
  result = ((double)mantissa) / pow(2.0, fracBits);
  return result;
  /*
  if (mantissa < condition){

    mantissa -= (1ull << (fracBits)); // take away the imply one
    return float_internal_to_double ((float_internal){0,-1, (unsigned long long)
  mantissa}); // FIX }else{

    mantissa -= (1ull << (fracBits+1));
    // printf("case2 mantissa_to_double %016llx \n",(unsigned long
  long)mantissa); return float_internal_to_double ((float_internal){0,0,
  (unsigned long long) mantissa}); // FIX
  }
  */
}

// mantissa with imply bit
//  double_to_mantissa works in a interval of [0.5,2.0)
unsigned __int128 double_to_mantissa(double x) {
  double_union_t ux;
  unsigned __int128 mantissa = 0;
  ux.d = x;
  // printf("exp %d %016llX ", fix.exp, fix.exp);
  int exp;
  exp = (ux.ull & 0x7FF0000000000000ull) >> 52;
  exp -= 1023; // subtract bias
  mantissa += (1ull << (fracBits + exp)) /*implybit*/;
  // int shift = x < 1.0 ?  FRACBITS - 1 : FRACBITS;
  mantissa += (fracBits + exp) > 52
                  ? (ux.ull & 0x000FFFFFFFFFFFFFull) << ((fracBits + exp) - 52)
                  : (ux.ull & 0x000FFFFFFFFFFFFFull) >>
                        (52 - (fracBits + exp)); // without assumed 1

  // printf("x %10.10f u %016llx exp %d %016x \n",x ,ux.ll,exp,exp);
  // fix = exp < 0 ? fix + (1ull << (FRACBITS-1)) : fix + (1ull << (FRACBITS));
  // // add imply 1  printf("fix %016llX \n",fix);
  return mantissa;
}

unsigned __int128 double_to_mantissa_withBit(double x, long long bits) {
  double_union_t ux;
  unsigned __int128 mantissa = 0;
  ux.d = x;
  // printf("exp %d %016llX ", fix.exp, fix.exp);
  int exp;
  exp = (ux.ull & 0x7FF0000000000000ull) >> 52;
  exp -= 1023; // subtract bias
  mantissa += (1ull << (bits + exp)) /*implybit*/;
  // int shift = x < 1.0 ?  FRACBITS - 1 : FRACBITS;
  mantissa += (bits + exp) > 52
                  ? (ux.ull & 0x000FFFFFFFFFFFFFull) << ((bits + exp) - 52)
                  : (ux.ull & 0x000FFFFFFFFFFFFFull) >>
                        (52 - (bits + exp)); // without assumed 1

  // printf("mantissa %016llx x %10.10f u %016llx exp %d %016x \n",mantissa,x
  // ,ux.ull,exp,exp);  fix = exp < 0 ? fix + (1ull << (FRACBITS-1)) : fix + (1ull
  // << (FRACBITS)); // add imply 1  printf("fix %016llX \n",fix);
  return mantissa;
}

unsigned long long float_to_fixed23(float x) {
  float_union_t ux;
  ux.f = x;
  int exp = (ux.i & 0x7F800000) >> 23;
  exp -= 127;
  unsigned long long significand = (1ull << (23 + exp)) /*implybit*/;
  significand +=
      exp > 0 ? (ux.i & 0x007FFFFF) << exp : (ux.i & 0x007FFFFF) >> (-exp);
  assert(!(exp < 0) || significand < (1ull << 23));
  return significand;
}

// Convert a float to fixed point precision with precision of 23 bits
unsigned long long float_to_fixed23_given_exp(float x, int unbiased_exp) {
  float_union_t ux;
  ux.f = x;
  int exp = (ux.i & 0x7F800000) >> 23;
  exp -= 127;
  assert(exp == unbiased_exp);
  unsigned long long significand = (1ull << 23) /*implybit*/;
  significand += ux.i & 0x007FFFFF;
  return significand;
}

// Convert a float to fixed point precision with specified precision bits
unsigned long long float_to_fixed23_withprec(float x, unsigned int prec) {
  float_union_t ux;
  unsigned long long significand;
  ux.f = x;
  int exp = (ux.i & 0x7F800000) >> 23;
  exp -= 127; // subtract bias
  unsigned long long implied_one =
      exp < 0 ? (1ull << 23) >> (-exp) : (1ull << 23) << exp;
  unsigned long long fraction = exp < 0
                                    ? (ux.i & 0x007FFFFF) >> (prec - 23 - exp)
                                    : (ux.i & 0x007FFFFF) << (prec - 23 + exp);
  // printf("implied_one: %016llx\n", implied_one);
  significand = fraction | implied_one;
  return significand;
}

int count_leading_zeros(unsigned long long x) {
  int n = 0;
  while (x < (1ull << 63)) {
    n++;
    x <<= 1;
  }
  return n;
}

float fixed23_to_float_given_exp(unsigned long long x, int unbiased_exp) {
  float result;
  result = ((float)x) / pow(2.0, 23 - unbiased_exp);
  return result;
}

// Convert internal float struct representation to IEEE float.
float float_internal_to_float(float_internal fix) {
  float x;
  if (fix.exp >= -126) {                 // normalized
    x = fix.significand * pow(2.0, -23); // fraction
    x += 1;                              // assumed 1
    x *= pow(2.0, fix.exp);              // exponent
  } else {                               // denormalized
    printf("float_internal_to_float:  denormalized x not implemented\n");
    assert(0);
  }
  return x;
}

double float_internal_to_double_mast(float_internal fix, int masBits) {
  double x;
  if (fix.exp >= -1022) {                         // normalized
    x = fix.significand * pow(2.0, -1 * masBits); // fraction
    x += 1;                                       // assumed 1
    x *= pow(2.0, fix.exp);                       // exponent
  } else {                                        // denormalized
    printf("float_internal_to_double:  denormalized x not implemented\n");
    assert(0);
  }
  return x;
}

// round to zero
float_internal double_to_float_internal(double x) {
  float_internal fix;
  double_union_t ux;
  ux.d = x;
  fix.sign = (x >= 0 ? 0 : 1);
  fix.exp = ((ux.ull & 0x7FF0000000000000ull) >> (64 - 12)) - 1023;
  // printf("exp %d %016llX ", fix.exp, fix.exp);
  fix.significand =
      (ux.ull & 0x000FFFFFFFFFFFFFull) >> (52 - FRACBITS); // without assumed 1
  // printf("ull  %016llX  signif %016llX \n",ux.ll , fix.significand);
  return fix;
}

// Convert IEEE float to internal float struct representation.
float_internal float_to_float_internal(float x) {
  float_internal fix;
  float_union_t ux;
  ux.f = x;
  fix.sign = (x >= 0 ? 0 : 1);
  fix.exp = ((ux.i & 0x7F800000) >> 23) - 127; // subtract IEEE bias
  fix.significand = ux.i & 0x007FFFFF;         // without assumed 1
  return fix;
}

double mypow2(int expo) {
  double res;
  // Compute 2**expo being careful about pow() problems with denorms
  if (expo >= -1022) {
    res = pow(2.0, expo);
  } else {
    res = pow(2.0, expo + 53) * pow(2.0, -53);
  }
  return res;
}

double ulp_d(double y) {
  long long expobits;
  int expo;
  union {
    double d;
    unsigned long long ull;
  } u;
  u.d = y;
  // printf("u.ull    =%08llX\n", u.ull);

  expobits = u.ull & 0x7FF0000000000000ull;
  // printf("expobits=%08X\n", expobits);
  expo = expobits >> (64 - 12);
  // printf("expo    =%08X=%d\n", expo, expo);
  expo -= 1023;
  // printf("expo    =%08X=%d\n", expo, expo);
  if (expo >= -1022) {
    return mypow2(expo - FRACBITS);
  } else {
    printf("ulp_d denormalized\n");
    assert(0);
    return mypow2(-149);
  }
}

float ulp_s(float y) {
  unsigned int expobits;
  int expo;
  union {
    float f;
    int i;
  } u;
  u.f = y;
  // printf("u.i    =%08X\n", u.i);
  expobits = u.i & 0x7F800000;
  // printf("expobits=%08X\n", expobits);
  expo = expobits >> (32 - 9);
  // printf("expo    =%08X=%d\n", expo, expo);
  expo -= 127;
  // printf("expo    =%08X=%d\n", expo, expo);
  if (expo >= -126) {
    return mypow2(expo - 23);
  } else {
    return mypow2(-149);
  }
}
