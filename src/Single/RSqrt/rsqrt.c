
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

#include "float_defs.h"
#include "float_utils.h"
#include "rsqrt.h"

// The T-Table generation uses round values of the function of interest
// This is the `rsqrtx` variable.
void makeTTable(__uint64_t* table) {
  int i;
  for (i=0;i<RSQRT_TABLE_SIZE/2;i++) {
    double num =  1 + oneOverSize * i;

    double rsqrtx = 1/sqrt(num);
    if (i==0) {
      table[i] =  (1ull<<(fracBitsTable+1)) << (fracBits-(fracBitsTable+1));
      table[256] = (__int64_t)(floor(1/sqrt(2) * (1ull<<fracBitsTable) * 2)) << (fracBits-(fracBitsTable+1));
    } else {
      table[i] = (__int64_t)(floor(rsqrtx * (1ull<<fracBitsTable) * 2)) << (fracBits-(fracBitsTable+1));
      num =  2 + oneOverSize * 2 * i;
      rsqrtx = 1/sqrt(num);
      table[i+256] = (__int64_t)(floor(rsqrtx * (1ull<<fracBitsTable) * 2)) << (fracBits-(fracBitsTable+1));
    }
  }
}

// The S-Table generation uses adjacent values in the T-Table.
// This assumes `makeTTable` is called first.
void makeSTable(__uint64_t* table, __uint64_t* T_Table) {

  table[0] = ((T_Table[0])
    - (T_Table[1]));
  for (int i=1;i<RSQRT_TABLE_SIZE-1;i++){
    table[i] = (T_Table[i]-T_Table[i+1]);
  }
  table[RSQRT_TABLE_SIZE-1] =
     ((T_Table[(RSQRT_TABLE_SIZE-1)])
      - ((1ull<<(fracBits-1))));
}

// This is a legacy function that calculates an approximation
// with just the T-Table and S-Table.
// This is used for the generation of the A-Table.
unsigned __int128 approxRSqrtLinear(unsigned __int128 x, __uint64_t* T_Table, __uint64_t* S_Table) {
    __uint64_t leadingBit = ((1ull << (fracBits+1)) & x) == 0 ? 0 : 1;
    __uint64_t newLowu = leadingBit == 0 ? lowu : lowu + 1;
    __uint64_t idx = ((1ull<<highf)-1) & (((__uint64_t) x) >> newLowu);
    idx += leadingBit << highf;

    signed __int128 lowBits = ((1ull<<newLowu)-1) & (x);

    // the first exponent in the table is 0 (for value 1.0) and -1 for the rest (1/2 < value < 1)
    bool idx_isOne = idx == 0 ;
    __uint64_t tableSignificand = (T_Table[idx]);

    signed __int128 S_Tableslope = (S_Table[idx]); //less than 20 bits

    S_Tableslope *= lowBits;

    signed __int128 l ;
    l = tableSignificand;// 50
    l <<= newLowu;


    l -= S_Tableslope;
    l >>= newLowu;

    return l;
}

// The Q-Table is generated via a table search method.
// For each subinterval, we take the worst-case correction value
// of `approxRecipLinear`.
void makeQTable (__uint64_t* table, __uint64_t* T_Table, __uint64_t* S_Table){
  // input x
  // oneOverSize = 1 / 256
  int count = 0;
  // first half range [1,2)
  for(int id = 0 ; id < RSQRT_TABLE_SIZE/2 ; id ++) {
    for (int alpha = 0 ; alpha < pow(2,qTableBit); alpha += 1){
      double x1 = 1
	+ ((id * pow(2,qTableBit)) + alpha)
	* (oneOverSize * (1/pow(2,qTableBit))); //1.000244140625
      double x2 = 1
	+ ((id * pow(2,qTableBit)) + (alpha + 1))
	* (oneOverSize * (1/pow(2,qTableBit))); //1.00048828125
      double step = pow (2, -23);
      __int64_t result = 0;

      for(double x = x1 ; x <= x2 ; x += step){

	double rsqrtx = 1/sqrt(x);
	unsigned __int128 tmp = double_to_mantissa (x);
	unsigned __int128 l_estimate = approxRSqrtLinear(tmp, T_Table, S_Table);
	double l = mantissa_to_double(l_estimate);

	__int64_t sub = (__int64_t)(floor((l - rsqrtx) * pow(2,23)));
	__int64_t tmpR = result;
	result = x == x1?/*init result*/ sub : (result > sub ? sub : result);
	if (result < 0) result = 0;

	if(tmpR != result){
	  //printf("x %16.16f l %16.16f rsqrtx %16.16f before floor %10.10f sub %lld result %lld %016llx\n",x,l,rsqrtx,(l - rsqrtx) * pow(2,23),sub,result,result);
	   }
	}// end of x

      table [count] = result;
      count ++;
    }// end of alpha
  } // end of id
  // second half range [2,4)
  for(int id = 0 ; id < RSQRT_TABLE_SIZE/2 ; id ++) {
    for (int alpha = 0 ; alpha < pow(2,qTableBit); alpha += 1){
      double x1 = 2
	+ ((id * pow(2,qTableBit)) + alpha)
	* (oneOverSize * 2 * (1/pow(2,qTableBit))); //1.000244140625
      double x2 = 2
	+ ((id * pow(2,qTableBit)) + (alpha + 1))
	* (oneOverSize * 2 * (1/pow(2,qTableBit))); //1.00048828125
      double step = pow (2, -23); //check
      __int64_t result = 0;

      for(double x = x1 ; x <= x2 ; x += step){

	double rsqrtx = 1/sqrt(x);
	unsigned __int128 tmp = double_to_mantissa (x);
	unsigned __int128 l_estimate = approxRSqrtLinear (tmp, T_Table, S_Table);
	double l = mantissa_to_double(l_estimate);

	__int64_t sub = (__int64_t)fabs(floor (((l - rsqrtx) * pow(2,23))));
	__int64_t tmpR = result;
	result = x == x1?/*init result*/ sub : (result > sub ? sub : result);
	if (result < 0) result = 0;

	if(tmpR != result){
	  //printf("x %16.16f l %16.16f rsqrtx %16.16f before floor %10.10f sub %lld result %lld %016llx\n",x,l,rsqrtx,(l - rsqrtx) * pow(2,23),sub,result,result);
	   }
	}// end of x

      table [count] = result;
      count ++;
    }// end of alpha
  } // end of id
  /*
  count = 0;
  printf("make q table \n");
  for(int id = 0 ; id < SQRT_TABLE_SIZE ; id ++) {
    for (int alpha = 0 ; alpha < pow(2,qTableBit); alpha += 1){
      printf ("(%i , %016llx)\n",alpha,table[count]);

      count ++;
    }
    printf("\n");
  }
  exit (0);
  */
}

// The A-Table is a subset of the Q-Table, and the entirety of the
// Q-Table can be restored via the bit-complement trick described in the paper.
void makeATable(__uint8_t *table, __uint64_t *qTable) {
  for (int i = 0; i < RSQRT_TABLE_SIZE; i++) {
    int qTableIdx = ((1ull << qTableBit) / 2) + (i * (1ull << qTableBit));
    table[i] = qTable[qTableIdx];
  }
}

// The main three-table method, which assumes that the T, S, and A-Table
// are generated.
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
  __uint32_t y = t - L - Q;

  return y;
}

double calculate_ytrue(unsigned __int128 x){
  return 1/sqrt(mantissa_to_double(x));
}

double calculate_ytest(unsigned __int128 x,__uint64_t* T_Table, __uint64_t* S_Table, __uint8_t* A_Table){
  return mantissa_to_double(three_table_procedure(x,T_Table,S_Table,A_Table));
}

int calculate_id(unsigned __int128 x){
  __uint64_t leadingBit = ((1ull << (fracBits+1)) & x) == 0 ? 0 : 1;
  __uint64_t newLowu = leadingBit == 0 ? lowu : lowu + 1;
  __uint64_t idx = ((1ull<<highf)-1) & (((__uint64_t) x) >> newLowu);
  idx += leadingBit << highf;
  return (int)idx;
}
