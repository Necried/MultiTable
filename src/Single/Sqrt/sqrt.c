//
//  main.c
//  DeepFractions
//
//  Created by Christopher Kumar Anand on 2021-01-27.
//


#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "float_defs.h"
#include "float_utils.h"
#include "sqrt.h"

// RECIP_TABLE_SIZE is (256)
// table makePosRecipTable with imply bits
void makeTTable(__uint64_t *table) {
  int i;
  for (i = 0; i < SQRT_TABLE_SIZE / 2; i++) {
    double num = 1 + oneOverSize * i;

    double sqrtx = sqrt(num);
    if (i == 0) {
      table[i] = (1ull << (fracBitsTable + 1))
                 << (fracBits - (fracBitsTable + 1));
      table[256] = (long long)(floor(sqrt(2) * (1ull << fracBitsTable) * 2))
                   << (fracBits - (fracBitsTable + 1));
    } else {
      // table[i] = (ceil((inv - 0.5 /*implied bit*/) * (1ll<<fracBitsTable) *
      // 2)) ;
      table[i] = (long long)(floor(sqrtx * (1ull << fracBitsTable) * 2))
                 << (fracBits - (fracBitsTable + 1));
      num = 2 + oneOverSize * 2 * i;
      sqrtx = sqrt(num);
      table[i + 256] = (long long)(floor(sqrtx * (1ull << fracBitsTable) * 2))
                       << (fracBits - (fracBitsTable + 1));
    }
    // printf("%i %016x %f %10.16f\n",i,table[i],num,inv);
  }
  // printf("makePosRecipTable \n");
  for (i = 0; i < SQRT_TABLE_SIZE; i++) {
    // printf(" (%i, %016llx, %f)\n",i,table[i],mantissa_to_double(table[i]));
  }
}

// slope
void makeSTable(__uint64_t *table, __uint64_t *sqrtable) {

  table[0] = (sqrtable[1]) - (sqrtable[0]);
  for (int i = 1; i < SQRT_TABLE_SIZE - 1; i++) {
    table[i] = (sqrtable[i + 1] - sqrtable[i]);
  }
  table[SQRT_TABLE_SIZE - 1] =
      ((1ull << (fracBits + 1)) - ((sqrtable[(SQRT_TABLE_SIZE - 1)])));
}

/*  bit values
    <sign> (2^31) <exp> (2^30)..(2^23) <also implied 1 bit> (2^23) <fraction>
   (2^22)...(2^2)(2=2^1)(1=2^0) if exp=0, interval [1,2) (1/2=2^{-1})  ...
   (2^{-23}) (2^(-23)*2^(-23) = 2^(-46) place value for ulp after multiply)
    after multiply, the place values are
    t2 <carry>(2^0)(2^(-1))(2^(-2))...(2^(-23))...(2^(-46))
 */
unsigned __int128 approxSqrtLinear(unsigned __int128 x, __uint64_t *T_Table,
                                   __uint64_t *S_Table) {
  __uint64_t leadingBit = ((1ull << (fracBits + 1)) & x) == 0 ? 0 : 1;
  __uint64_t newLowu = leadingBit == 0 ? lowu : lowu + 1;
  __uint64_t idx =
      ((1ull << highf) - 1) & (((__uint64_t)x) >> newLowu);
  idx += leadingBit << highf;

  unsigned __int128 lowBits = ((1ull << newLowu) - 1) & (x); // 16 bits

  // the first exponent in the table is 0 (for value 1.0) and -1 for the rest
  // (1/2 < value < 1)
  bool idx_isOne = idx == 0;
  __uint64_t tableSignificand = 0;
  __uint64_t table = (T_Table[idx]);

  tableSignificand = table;

  __uint64_t S_Tableslope = (S_Table[idx]); // less than 20 bits

  unsigned __int128 l;
  l = tableSignificand; // 50
  l += (S_Tableslope * lowBits) >> newLowu;

  // printf("approxRecipLiner \n");
  // printf("x 0x%016llx %016llx idx %d lowBits 0x%016llx l 0x%016llx 0x%016llx
  // tableSignificand 0x%016llx tableSignificand(from approxRecip) 0x%016llx
  // S_Tableslope 0x%016llx \n",(__uint64_t)(x>>64),(unsigned long
  // long)x,idx,(__uint64_t)lowBits,(__uint64_t)
  // (l>>64),(__uint64_t) (l),(__uint64_t)
  // (tableSignificand),recipTable[idx],S_Tableslope);

  // printf("in function approxRecipLiner %10.16f %016llx
  // \n",mantissa_to_double(l),l);
  return l;
}

void makeQTable(__uint64_t *table, __uint64_t *T_Table,
                __uint64_t *S_Table) {
  // input x
  // oneOverSize = 1 / 256
  int count = 0;
  // first half range [1,2)
  for (int id = 0; id < SQRT_TABLE_SIZE / 2; id++) {
    for (int alpha = 0; alpha < pow(2, qTableBit); alpha += 1) {
      double x1 =
          1 + ((id * pow(2, qTableBit)) + alpha) *
                  (oneOverSize * (1 / pow(2, qTableBit))); // 1.000244140625
      double x2 =
          1 + ((id * pow(2, qTableBit)) + (alpha + 1)) *
                  (oneOverSize * (1 / pow(2, qTableBit))); // 1.00048828125
      double step = pow(2, -23);
      long long result = 0;

      for (double x = x1; x <= x2; x += step) {

        double sqrtx = sqrt(x);
        unsigned __int128 tmp = double_to_mantissa(x);
        unsigned __int128 l_estimate =
            approxSqrtLinear(tmp, T_Table, S_Table);
        double l = mantissa_to_double(l_estimate);

        long long sub = (long long)floor((sqrtx - l) * pow(2, 23));
        long long tmpR = result;
        result = x == x1 ? /*init result*/ sub : (result > sub ? sub : result);

        if (tmpR != result) {
          // printf("x %16.16f l_double %16.16f inv %16.16f before floor %10.10f
          // sub %d result %d %016x\n",x,l_double,inv,(l_double - inv) *
          // pow(2,23),sub,result,result);
        }
      } // end of x

      table[count] = result;
      count++;
    } // end of alpha
  }   // end of id
  // second half range [2,4)
  for (int id = 0; id < SQRT_TABLE_SIZE / 2; id++) {
    for (int alpha = 0; alpha < pow(2, qTableBit); alpha += 1) {
      double x1 =
          2 + ((id * pow(2, qTableBit)) + alpha) *
                  (oneOverSize * 2 * (1 / pow(2, qTableBit))); // 1.000244140625
      double x2 =
          2 + ((id * pow(2, qTableBit)) + (alpha + 1)) *
                  (oneOverSize * 2 * (1 / pow(2, qTableBit))); // 1.00048828125
      double step = pow(2, -23);                               // check
      long long result = 0;

      for (double x = x1; x <= x2; x += step) {

        double sqrtx = sqrt(x);
        unsigned __int128 tmp = double_to_mantissa(x);
        unsigned __int128 l_estimate =
            approxSqrtLinear(tmp, T_Table, S_Table);
        double l = mantissa_to_double(l_estimate);

        long long sub = (long long)floor(((sqrtx - l) * pow(2, 23)));
        long long tmpR = result;
        result = x == x1 ? /*init result*/ sub : (result > sub ? sub : result);

        if (tmpR != result) {
          // printf("x %16.16f l_double %16.16f inv %16.16f before floor %10.10f
          // sub %d result %d %016x\n",x,l_double,inv,(l_double - inv) *
          // pow(2,23),sub,result,result);
        }
      } // end of x

      table[count] = result;
      count++;
    } // end of alpha
  }   // end of id
  count = 0;
  // printf("make q table \n");
  for (int id = 0; id < SQRT_TABLE_SIZE; id++) {
    for (int alpha = 0; alpha < pow(2, qTableBit); alpha += 1) {
      // printf ("(%i , %016llx)\n",alpha,table[count]);

      count++;
    }
    // printf("\n");
  }
  // exit (0);
}

void makeATable(__uint8_t *table, __uint64_t *qTable) {
  for (int i = 0; i < SQRT_TABLE_SIZE; i++) {
    int qTableIdx = ((1ull << qTableBit) / 2) + (i * (1ull << qTableBit));
    table[i] = qTable[qTableIdx];
  }
}

unsigned __int128 three_table_procedure(unsigned __int128 x,
                                        __uint64_t *T_Table,
                                        __uint64_t *S_Table,
                                        __uint8_t *A_Table) {
  // Extract indices of x
  __uint32_t leadingBit = ((1ull << (fracBits + 1)) & x) == 0 ? 0 : 1;
  __uint32_t newLowu = leadingBit == 0 ? lowu : lowu + 1;
  __uint32_t x0__8 = (leadingBit << highf) + (((1ull << highf) - 1) &
                     (((__uint64_t)x) >> newLowu));
  assert(0 <= x0__8 && x0__8 < 512);

  __uint32_t x9__18 =
      ((1ull << qTableBit) - 1) & (((__uint64_t)x) >> (newLowu - qTableBit));
  __uint32_t x9__23 = ((1ull << newLowu) - 1) & (x);


  // Perform table lookups
  __uint32_t t = T_Table[x0__8];
  __uint32_t s = S_Table[x0__8];
  __uint32_t a = A_Table[x0__8];
  
  // Calculate z
  // Here, we first compute the complement of x[9..18] separately
  // The complement operator in C produces extraneous leading bits, so we need
  // to mask.
  __uint32_t x9__18_comp = ((1 << qTableBit) - 1) & (~x9__18);
  __uint32_t z = x9__18 * x9__18_comp;
  // Ensures that the result of the multiplication is at most 22 bits wide
  // assert(z < (1 << 22));

  // The base shift amount of Q is 18, but we take into account extra precision
  // introduced by fracBits. We subtract by 23 (the single precision
  // significand)
  __uint32_t Q = (z * a) >> (18 - (fracBits - 23));
  __uint32_t L = ((__uint64_t)s * x9__23) >> newLowu;
  // The difference from recip is we add the correction instead of subtract.
  __uint32_t y = t + L + Q;

  return y;
}

double calculate_ytrue(unsigned __int128 x) {
  return sqrt(mantissa_to_double(x));
}

double calculate_ytest(unsigned __int128 x,
                       __uint64_t *T_Table, __uint64_t *S_Table,
                       __uint8_t *A_Table) {
  return mantissa_to_double(
      three_table_procedure(x, T_Table, S_Table, A_Table));
}

int calculate_id(unsigned __int128 x) {
  __uint64_t leadingBit = ((1ull << (fracBits + 1)) & x) == 0 ? 0 : 1;
  __uint64_t newLowu = leadingBit == 0 ? lowu : lowu + 1;
  __uint64_t idx =
      ((1ull << highf) - 1) & (((__uint64_t)x) >> newLowu);
  idx += leadingBit << highf;
  return (int)idx;
}