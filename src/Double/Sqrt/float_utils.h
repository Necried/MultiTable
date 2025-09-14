#ifndef DOUBLE_RECIP_FLOAT_UTILS_H
#define DOUBLE_RECIP_FLOAT_UTILS_H

#define FORMAT_UINT128(x) (unsigned long long)(x>>64), (unsigned long long)(x)
#define FORMAT_UINT128_S(x) STR(x), (unsigned long long)(x>>64), (unsigned long long)(x)

extern double mantissa_to_double (unsigned __int128 mantissa);
extern unsigned __int128 double_to_mantissa (double x);

double ieee_to_double(unsigned long long x);

extern double mypow2 (int expo);

extern double ulp_d (double y);
extern float ulp_s (float y);

extern unsigned __int128 big_mult(unsigned __int128 a, unsigned __int128 b);
#endif

#define PRINT_FORMAT_UINT128   "%-20s: %016llx %016llx\n"
#define PRINT_FORMAT_UINT128_2 "%-20s: %016llx %016llx\n%-20s: %016llx %016llx\n"
#define PRINT_FORMAT_UINT128_3 "%-20s: %016llx %016llx\n%-20s: %016llx %016llx\n%-20s: %016llx %016llx\n"
#define PRINT_FORMAT_UINT128_4 "%-20s: %016llx %016llx\n%-20s: %016llx %016llx\n%-20s: %016llx %016llx\n%-20s: %016llx %016llx\n"

#define DEBUG_UINT128(x) DEBUG_PRINT((PRINT_FORMAT_UINT128, FORMAT_UINT128_S(x)))
