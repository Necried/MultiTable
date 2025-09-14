#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h> // Check available stack space
#include <time.h>

#include "dtof_rounding.h"
#include "float_defs.h"
#include "float_utils.h"

#include "rsqrt.h"
#include "rsqrt_dp.h"

#define test_count TESTCOUNT

double gen_random_double_1_4() {
  union {unsigned long long ull; double d; int i;} val;
  assert(RAND_MAX == 2147483647);

  // unsigned long long top31 = ((unsigned long long)rand() % (1<<23)-1) << (52-23);
  unsigned long long top31 = (unsigned long long)rand() << (52-31);
  unsigned long long bot21 = (unsigned long long)rand() % ((1<<21)-1);
  unsigned long long rand_binade = (unsigned long long)rand() % 2;
  unsigned long long leadingOne = ((1ull << 10) - 1 + rand_binade) << 52;

  val.ull = leadingOne + top31 + bot21;
  // printf("top31 0x%llx   bot21 0x%llx\nleadingOne 0x%llx   val.ull 0x%llx\n",top31,bot21,leadingOne,val.ull);
  // printf("Generated double: %10.10f\n", val.d);
  assert(1 <= val.d && val.d < 4);

  return val.d;
}

void test_rand_doubles(unsigned int seed, int ontests, __uint64_t *T_Table,
          __uint64_t *S_Table, __uint8_t *A_Table) {
  // Initialize historgram vars
  long long n_ulp_05      = 0;
  long long n_ulp_1       = 0;
  long long n_ulp_2       = 0;
  long long n_ulp_5       = 0;
  long long n_ulp_10      = 0;
  long long n_ulp_10up    = 0;
  long long n_ulp_n05     = 0;
  long long n_ulp_n1      = 0;
  long long n_ulp_n2      = 0;
  long long n_ulp_n5      = 0;
  long long n_ulp_n10     = 0;
  long long n_ulp_n10down = 0;
  long long ncorrnd       = 0;

  double ytrue, ytest, d, ulperr;
  double omaxulperr = 0, omaxx = 0;
  double ominulperr = 0, ominx = 0;

  double omaxieee, ominieee;

  union {unsigned long long ull; double d;}
       omaxytest, omaxytrue, ominytest, ominytrue;

  srand(seed);

  for (int i = 0; i < ontests; i++) {
    d = gen_random_double_1_4();
    unsigned __int128 x = double_to_mantissa(d);
    ytrue = calculate_ytrue_dp(x);
    ytest = calculate_ytest_dp_newton(x, T_Table, S_Table, A_Table);
    ulperr = (ytest - ytrue) / (ulp_d (ytrue));
    //update histogram
    if(ulperr<-10){
	 n_ulp_n10down++;
    }else if(ulperr<-5.0){
	 n_ulp_n10++;
    }else if(ulperr<-2.0){
	 n_ulp_n5++;
    }else if(ulperr<-1.0){
	 n_ulp_n2++;
    }else if(ulperr<-0.5){
	 n_ulp_n1++;
    }else if(ulperr<=0){
	 n_ulp_n05++;
	 // start counting positive
    }else if(ulperr<=0.5){
	 n_ulp_05 ++;
    }else if (ulperr<=1){
	 n_ulp_1 ++;
    }else if (ulperr<=2){
	 n_ulp_2 ++;
    }else if (ulperr<=5){
	 n_ulp_5 ++;
    }else if (ulperr<=10){
	 n_ulp_10 ++;
    }else{
	 n_ulp_10up ++;
    }

    if(fabs(ulperr)<=0.5){
	 ncorrnd++;
    }

    // Get overall max and min
    if(ulperr > omaxulperr){
	 omaxieee = d;
	 omaxulperr = ulperr;
	 omaxytest.d = ytest;
	 omaxytrue.d = ytrue;
	 omaxx = x;
    }
    if(ulperr < ominulperr){
	 ominieee = d;
	 ominulperr = ulperr;
	 ominytest.d = ytest;
	 ominytrue.d = ytrue;
	 ominx = x;
    }
  }

  printf (" == ieeemant ==       == x ==   == frac ==      ========== ytest =========   ========= ytrue ========     ulps\n");
  printf("Max \n");

  printf ("%10.16f %10.2e %08llX %10.2e %016llX %10.2e %016llX %10.2e\n",
	  omaxieee,
	  mantissa_to_double(omaxx), (unsigned long long)omaxx,
	  omaxytest.d, omaxytest.ull,
	  omaxytrue.d, omaxytrue.ull,
	  omaxulperr);

  printf("Min \n");

  printf ("%10.16f %10.2e %08llX %10.2e %016llX %10.2e %016llX %10.2e\n",
	  ominieee,
	  mantissa_to_double(ominx), (unsigned long long)ominx,
	  ominytest.d, ominytest.ull,
	  ominytrue.d, ominytrue.ull,
	  ominulperr);
  printf("\n");

  printf("\n");
  printf("Number of test %d\n",ontests);
  printf ("===neg to pos===     ===count===  =====absolute===== \n");
  printf ("[-inf,-10) %6.2f%%  %10lld\n", 100.0*n_ulp_n10down/ontests,n_ulp_n10down);
  printf ("  [-10,-5) %6.2f%%  %10lld\n", 100.0*n_ulp_n10/ontests,n_ulp_n10);
  printf ("   [-5,-2) %6.2f%%  %10lld\n", 100.0*n_ulp_n5/ontests,n_ulp_n5);
  printf ("   [-2,-1) %6.2f%%  %10lld\n", 100.0*n_ulp_n2/ontests,n_ulp_n2);
  printf (" [-1,-0.5) %6.2f%%  %10lld\n", 100.0*n_ulp_n1/ontests,n_ulp_n1);
  printf ("  (-0.5,0] %6.2f%%  %10lld\n", 100.0*n_ulp_n05/ontests,n_ulp_n05);
  printf ("   (0,0.5] %6.2f%%  %10lld     [0,0.5] %6.2f%%\n", 100.0*n_ulp_05/ontests
						 ,n_ulp_05
						 , 100.0*(n_ulp_05+n_ulp_n05)/ontests );
  printf ("   (0.5,1] %6.2f%%  %10lld     (0.5,1] %6.2f%%\n", 100.0*n_ulp_1/ontests
						 ,n_ulp_1
						 , 100.0*(n_ulp_1+n_ulp_n1)/ontests);
  printf ("     (1,2] %6.2f%%  %10lld       (1,2] %6.2f%%\n", 100.0*n_ulp_2/ontests
						 ,n_ulp_2
						 , 100.0*(n_ulp_2+n_ulp_n2)/ontests);
  printf ("     (2,5] %6.2f%%  %10lld       (2,5] %6.2f%%\n", 100.0*n_ulp_5/ontests
						 ,n_ulp_5
						 , 100.0*(n_ulp_5+n_ulp_n5)/ontests);
  printf ("    (5,10] %6.2f%%  %10lld      (5,10] %6.2f%%\n", 100.0*n_ulp_10/ontests
						 ,n_ulp_10
						 , 100.0*(n_ulp_10+n_ulp_n10)/ontests);
  printf ("  (10,inf] %6.2f%%  %10lld    (10,inf] %6.2f%%\n", 100.0*n_ulp_10up/ontests
						 ,n_ulp_10up
						 , 100.0*(n_ulp_10up+n_ulp_n10down)/ontests);
  printf ("\n");
}

int main() {
  __uint64_t T_Table[RSQRT_TABLE_SIZE];
  __uint64_t S_Table[RSQRT_TABLE_SIZE];
  __uint64_t *qTable;
  __uint8_t A_Table[RSQRT_TABLE_SIZE];

  // Allocate heap space for qTable
  printf("QTableSize: %d\n", QTableSize);
  qTable = malloc(sizeof(__uint64_t) * QTableSize);

  makeTTable(T_Table);
  makeSTable(S_Table, T_Table);
  makeQTable(qTable, T_Table, S_Table);
  makeATable(A_Table, qTable);
  
  // Seed random test
  time_t seed = time(NULL);
  printf("Seed: %ld\n", seed);
  printf("Test count: %d\n", test_count);

  test_rand_doubles(seed, test_count, T_Table, S_Table, A_Table);

  printf("Testing sign 1 bits, exp 8 bits, significand %d bits \n", fracBits);
  printf("Testing rsqrt table : %d bits table, size %d \n", fracBitsTable,
         RSQRT_TABLE_SIZE);
  printf("Testing input(highf) %d bits \n", highf);
  printf("Testing lowu %d bits \n", lowu);

  free(qTable);
  return 0;
}