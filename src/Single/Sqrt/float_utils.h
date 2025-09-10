#include "float_defs.h"

extern float float_internal_to_float (float_internal fix);
extern float_internal float_to_float_internal (float x);

extern double float_internal_to_double (float_internal fix);
extern double float_internal_to_double_mast (float_internal fix, int masBits);
extern float_internal double_to_float_internal (double x);

extern unsigned long long float_to_fixed23 (float x);
extern unsigned long long float_to_fixed23_given_exp(float x, int unbiased_exp);
extern unsigned long long float_to_fixed23_withprec(float x, unsigned int prec);
extern float fixed23_to_float_given_exp(unsigned long long int x, int unbiased_exp);

extern unsigned __int128 double_to_mantissa_withBit (double x,long long bits);

extern double mantissa_to_double (unsigned __int128 mantissa);
extern unsigned __int128 double_to_mantissa (double x);

extern double mypow2 (int expo);

extern double ulp_d (double y);
extern float ulp_s (float y);
