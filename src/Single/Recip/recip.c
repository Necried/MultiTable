#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "float_defs.h"
#include "float_utils.h"
#include "recip.h"

// The T-Table generation uses round values of the function of interest
// This is the `inv` variable.
void makeTTable(__uint64_t *table) {
  int i;
  for (i = 0; i < RECIP_TABLE_SIZE; i++) {
    double num = 1 + oneOverSize * i;
    double inv = 1 / num;
    if (i == 0) {
      table[i] = (1ull << (fracBitsTable + 1))
                 << (fracBits - (fracBitsTable + 1));
    } else {
      table[i] = (__uint64_t)(ceil((inv) * (1ull << fracBitsTable) * 2))
                 << (fracBits - (fracBitsTable + 1));
    }
  }
  
  return;
}

// The S-Table generation uses adjacent values in the T-Table.
// This assumes `makeTTable` is called first.
void makeSTable(__uint64_t *table, __uint64_t *reciptable) {

  table[0] = (reciptable[0]) - (reciptable[1]);
  for (int i = 1; i < RECIP_TABLE_SIZE - 1; i++) {
    table[i] = (reciptable[i] - reciptable[i + 1]);
  }
  table[RECIP_TABLE_SIZE - 1] =
      ((reciptable[(RECIP_TABLE_SIZE - 1)]) - (1ull << (fracBits - 1)));

  return;
}

// This is a legacy function that calculates an approximation
// with just the T-Table and S-Table.
// This is used for the generation of the A-Table.
unsigned __int128 approxRecipLinear(unsigned __int128 x, __uint64_t *recipTable,
                                   __uint64_t *linearTable) {
  __uint64_t idx =
      ((1ull << highf) - 1) & (((__uint64_t)x) >> lowu);
  unsigned __int128 lowBits = ((1ull << lowu) - 1) & (x); // 16 bits

  // the first exponent in the table is 0 (for value 1.0) and -1 for the rest
  // (1/2 < value < 1)
  bool idx_isOne = idx == 0;
  __uint64_t tableSignificand = recipTable[idx];

  __uint64_t lineartableslope = (linearTable[idx]); // less than 20 bits

  unsigned __int128 l;
  l = tableSignificand; // 50
  l -= (lineartableslope * lowBits) >> lowu;

  return l;
}


// The Q-Table is generated via a table search method.
// For each subinterval, we take the worst-case correction value
// of `approxRecipLinear`.
void makeQTable(__uint64_t *table, __uint64_t *recipTable,
                __uint64_t *linearTable) {
  // input x
  // oneOverSize = 1 / 256
  int count = 0;
  for (int id = 0; id < RECIP_TABLE_SIZE; id++) {
    for (int alpha = 0; alpha < pow(2, qTableBit); alpha += 1) {
      double x1 =
          1 + ((id * pow(2, qTableBit)) + alpha) *
                  (oneOverSize * (1 / pow(2, qTableBit))); // 1.000244140625
      double x2 =
          1 + ((id * pow(2, qTableBit)) + (alpha + 1)) *
                  (oneOverSize * (1 / pow(2, qTableBit))); // 1.00048828125
      double step = pow(2, -23);
      __uint64_t result = 0;

      for (double x = x1; x <= x2; x += step) {

        double inv = 1 / x;
        unsigned __int128 tmp = double_to_mantissa(x);
        unsigned __int128 l = approxRecipLinear(tmp, recipTable, linearTable);
        double l_double = mantissa_to_double(l);

        __uint64_t sub = (__uint64_t)floor((l_double - inv) * pow(2, 23));
        __uint64_t tmpR = result;
        result = x == x1 ? /*init result*/ sub : (result > sub ? sub : result);

      } // end of x

      table[count] = result;
      count++;
    } // end of alpha
  }   // end of id
}

// The A-Table is a subset of the Q-Table, and the entirety of the
// Q-Table can be restored via the bit-complement trick described in the paper.
void makeATable(__uint8_t *table, __uint64_t *qTable) {
  for (int i = 0; i < RECIP_TABLE_SIZE; i++) {
    int qTableIdx = ((1ull << qTableBit) / 2) + (i * (1ull << qTableBit));
    table[i] = qTable[qTableIdx];
  }
}

// The main three-table method, which assumes that the T, S, and A-Table
// are generated.
unsigned __int128 three_table_procedure(unsigned __int128 x, __uint64_t *T_Table,
                             __uint64_t *S_Table, __uint8_t *A_Table) {
  // Extract indices of x
  __uint32_t x1__8  = ((1ull << highf) - 1) & (((__uint64_t)x) >> lowu);;
  __uint32_t x9__20 = ((1ull << qTableBit) - 1) & (((__uint64_t)x) >> (lowu - qTableBit));
  __uint32_t x9__23 =  ((1ull << lowu) - 1) & (x);

  // Perform table lookups
  __uint32_t t = T_Table[x1__8];
  __uint32_t s = S_Table[x1__8];
  __uint32_t a = A_Table[x1__8];
  
  // Calculate z
  // Here, we first compute the complement of x[9..20] separately
  // The complement operator in C produces extraneous leading bits, so we need 
  // to mask.
  __uint32_t x9__20_comp = ((1 << qTableBit) - 1) & (~x9__20);
  __uint32_t z = x9__20 * x9__20_comp;
  // Ensures that the result of the multiplication is at most 22 bits wide
  assert(z < (1 << 22));
  
  // The base shift amount of Q is 22, but we take into account extra precision
  // introduced by fracBits. We subtract by 23 (the single precision significand)
  __uint32_t Q = fracBits - 23 >= 0 ?
    (z * a) >> (22 - (fracBits - 23)) :
    (z * a) << ((fracBits - 23) - 22);
  __uint32_t L = ((__uint64_t)s * x9__23) >> lowu;
  __uint32_t y = t - L - Q;

  return y;
}

double calculate_ytrue(unsigned __int128 x) {
  return 1 / (mantissa_to_double(x));
}

double calculate_ytest(unsigned __int128 x, __uint64_t *recipTable,
                       __uint64_t *linearTable, __uint8_t *aTable) {
  // 3-table implementation
  double three_table_result =
      mantissa_to_double(three_table_procedure(x, recipTable, linearTable, aTable));
  return three_table_result;
}

int calculate_id(unsigned __int128 x) {

  return ((1ull << highf) - 1) & (((__uint64_t)x) >> lowu);
}