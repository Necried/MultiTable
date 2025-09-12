#include <stdlib.h>
double calculate_ytest_dp(unsigned __int128 b, __uint64_t *recipTable,
                       __uint64_t *linearTable, __uint8_t *aTable);
double calculate_ytrue_dp(unsigned __int128 x);
double calculate_ytest_dp_newton(unsigned __int128 a, __uint64_t *recipTable,
                       __uint64_t *linearTable, __uint8_t *aTable);