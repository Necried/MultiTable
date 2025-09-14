#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h> // Check available stack space
#include <unistd.h>

#include "dtof_rounding.h"
#include "float_defs.h"
#include "float_utils.h"

#include "recip.h"

// statistic
char exc_bits[6], oexc_bits[6]; // *<>nz
void clear_exc_bits() { strcpy(exc_bits, "....."); }

void clear_oexc_bits() { strcpy(oexc_bits, "....."); }

void set_exc_ninf() {
  exc_bits[0] = '<';
  oexc_bits[0] = '<';
}

void set_exc_pinf() {
  exc_bits[1] = '>';
  oexc_bits[1] = '>';
}

void set_exc_nan() {
  exc_bits[2] = 'n';
  oexc_bits[2] = 'n';
}

void set_exc_zerr() {
  exc_bits[3] = 'z';
  oexc_bits[3] = 'z';
}

void set_exc_bad() {
  exc_bits[4] = '*';
  oexc_bits[4] = '*';
}

void test(int exp, int fromtable, int tableend, __uint64_t *T_Table,
          __uint64_t *S_Table, __uint8_t *A_Table, FILE *fptr) {
  // define variant
  int id, maxid, minid;
  __int64_t ncorrnd, omaxncorrnd, ominncorrnd, ontests, ntests, oontests,
      n_ulp_05, n_ulp_1, n_ulp_2, n_ulp_5, n_ulp_10, n_ulp_10up, n_ulp_n05,
      n_ulp_n1, n_ulp_n2, n_ulp_n5, n_ulp_n10, n_ulp_n10down, on_ulp_05,
      on_ulp_1, on_ulp_2, on_ulp_5, on_ulp_10, on_ulp_10up, on_ulp_n05,
      on_ulp_n1, on_ulp_n2, on_ulp_n5, on_ulp_n10, on_ulp_n10down;
  ncorrnd = 0;
  omaxncorrnd = 0;
  ominncorrnd = 0;
  ontests = 0;
  ntests = 0;
  oontests = 0;
  n_ulp_05 = 0;
  n_ulp_1 = 0;
  n_ulp_2 = 0;
  n_ulp_5 = 0;
  n_ulp_10 = 0;
  n_ulp_10up = 0;
  n_ulp_n05 = 0;
  n_ulp_n1 = 0;
  n_ulp_n2 = 0;
  n_ulp_n5 = 0;
  n_ulp_n10 = 0;
  n_ulp_n10down = 0;
  on_ulp_05 = 0;
  on_ulp_1 = 0;
  on_ulp_2 = 0;
  on_ulp_5 = 0;
  on_ulp_10 = 0;
  on_ulp_10up = 0;
  on_ulp_n05 = 0;
  on_ulp_n1 = 0;
  on_ulp_n2 = 0;
  on_ulp_n5 = 0;
  on_ulp_n10 = 0;
  on_ulp_n10down = 0;
  id = 0;
  maxid = 0;
  minid = 0;

  union {
    __uint64_t ull;
    double d;
    int i;
  } ytrue, ytest, maxytrue, maxytest, omaxytrue, omaxytest, ominytrue,
      ominytest;

  ytrue.d = 0;
  ytest.d = 0;
  maxytrue.d = 0;
  maxytest.d = 0;

  omaxytest.d = 0;
  ominytest.d = 0;
  omaxytrue.d = 0;
  ominytrue.d = 0;

  double maxulperr = 0;

  double omaxulperr = 0;
  double ominulperr = 0;

  unsigned __int128 maxx = 0;
  unsigned __int128 minx = 0;

  unsigned __int128 omaxx = 0;
  unsigned __int128 ominx = 0;

  // print header
  fprintf(fptr,"go through every significand from table %d to table %d \n", fromtable,
         tableend);
  fprintf(fptr,"go through exponent 0 to exponent %d\n", exp);
  fprintf(fptr,"increment 1 bit every line\n");
  fprintf(fptr,"\n");

  // go through exponent from 0 to exp
  for (int j = 0; j <= exp; j++) {
    // initial histogram
    n_ulp_05 = 0;
    n_ulp_1 = 0;
    n_ulp_2 = 0;
    n_ulp_5 = 0;
    n_ulp_10 = 0;
    n_ulp_10up = 0;
    n_ulp_n05 = 0;
    n_ulp_n1 = 0;
    n_ulp_n2 = 0;
    n_ulp_n5 = 0;
    n_ulp_n10 = 0;
    n_ulp_n10down = 0;
    ontests = 0;
    // go throught table fromtable to tableend
#ifdef TABLE
    fprintf(fptr,"id  sgn  exp <>nz*  == x == == frac == ========== ytest =========  "
           " ========= ytrue ========     ulps      iulps    pcr\n");
#endif
    for (int k = fromtable; k <= tableend; k++) {
      // initial for every table
      ncorrnd = 0;
      ntests = 0;

      maxulperr = 0;

      // maxid = 0;
      // go through every significand in a table interval
      __uint64_t at;
      __uint64_t i;

      for (i = k * (1l << lowu); i < (k + 1) * (1l << lowu);
           i += pow(2, fracBits - 23)) {
        at = k * (1l << lowu);

        // test
        unsigned __int128 x = (1ull << (fracBits)) + i;

        ytrue.d = calculate_ytrue(x);
        ytest.d = calculate_ytest(x, T_Table, S_Table, A_Table);
        double ulperr;
#ifdef SINGLE

#ifdef RU
        ytest.d = ((double)dtof_pinf(ytest.d));
        ytrue.d = ((double)dtof_pinf(ytrue.d));
#endif
#ifdef RD
        ytest.d = ((double)dtof_ninf(ytest.d));
        ytrue.d = ((double)dtof_ninf(ytrue.d));
#endif
#ifdef RZ
        ytest.d = ((double)dtof_zero(ytest.d));
        ytrue.d = ((double)dtof_zero(ytrue.d));
#endif
#ifdef RN
        ytest.d = ((double)dtof_nearest(ytest.d));
        ytrue.d = ((double)dtof_nearest(ytrue.d));
#endif
        ulperr = (ytest.d - ytrue.d) / (double)(ulp_s((float)ytrue.d));
#else
        ulperr = (ytest.d - ytrue.d) / (ulp_d(ytrue.d));
#endif

        id = calculate_id(x);

        // update histogram
        if (ulperr < -10) {
          n_ulp_n10down++;
          on_ulp_n10down++;
        } else if (ulperr < -5.0) {
          n_ulp_n10++;
          on_ulp_n10++;
        } else if (ulperr < -2.0) {
          n_ulp_n5++;
          on_ulp_n5++;
        } else if (ulperr < -1.0) {
          n_ulp_n2++;
          on_ulp_n2++;
        } else if (ulperr < -0.5) {
          n_ulp_n1++;
          on_ulp_n1++;
        } else if (ulperr <= 0) {
          n_ulp_n05++;
          on_ulp_n05++;
          // start counting positive
        } else if (ulperr <= 0.5) {
          n_ulp_05++;
          on_ulp_05++;
        } else if (ulperr <= 1) {
          n_ulp_1++;
          on_ulp_1++;
        } else if (ulperr <= 2) {
          n_ulp_2++;
          on_ulp_2++;
        } else if (ulperr <= 5) {
          n_ulp_5++;
          on_ulp_5++;
        } else if (ulperr <= 10) {
          n_ulp_10++;
          on_ulp_10++;
        } else {
          n_ulp_10up++;
          on_ulp_10up++;
        }

        if (fabs(ulperr) <= 0.5) {
          ncorrnd++;
        }
        ntests++;
        ontests++;
        oontests++;

        // get max and min in every table
        if (i == at) {
          maxulperr = fabs(ulperr);
          maxytest = ytest;
          maxytrue = ytrue;
          maxx = x;
        }

        if (fabs(ulperr) > fabs(maxulperr)) {
          maxulperr = fabs(ulperr);
          maxytest = ytest;
          maxytrue = ytrue;
          maxx = x;
        }

        // Get overall max and min
        if (ulperr > omaxulperr) {
          omaxulperr = ulperr;
          omaxytest = ytest;
          omaxytrue = ytrue;
          omaxx = x;
          maxid = id;
        }
        if (ulperr < ominulperr) {
          ominulperr = ulperr;
          ominytest = ytest;
          ominytrue = ytrue;
          ominx = x;
          minid = id;
        }

      }

        // print statistic for each table

#ifdef TABLE
      __int64_t maxiulperr = maxytest.i - maxytrue.i;
      fprintf(fptr,"%3d %2d %5d %5s %10.2e %010lX %10.2e %016lX %10.2e %016lX "
             "%10.2e %10.2e %6.2f\n",
             id, 0, 0, "     ",
             mantissa_to_double(maxx), (__int64_t)maxx, maxytest.d,
             maxytest.ull, maxytrue.d, maxytrue.ull, maxulperr,
             (double)maxiulperr, 100.0 * ncorrnd / ntests);
#endif

      if (maxid == id) {
        omaxncorrnd = ncorrnd;
      }
      if (minid == id) {
        ominncorrnd = ncorrnd;
      }

      // end k
    }

    // print statistic for each exp
    fprintf(fptr,"\n");
    fprintf(fptr,"overall of exponent %d\n", j);
    fprintf(fptr,"id  sgn  exp <>NZ*  == x == == frac == ========== ytest =========  "
           " ========= ytrue ========     ulps      iulps    pcr\n");
    fprintf(fptr,"Max \n");

    __int64_t omaxiulperr = omaxytest.i - omaxytrue.i;
    fprintf(fptr,"%3d %2d %5d %5s %10.2e %08lX %10.2e %016lX %10.2e %016lX %10.2e "
           "%10.2e %6.2f\n",
           maxid, 0, 0, "     ",
           mantissa_to_double(omaxx), (__int64_t)omaxx, omaxytest.d,
           omaxytest.ull, omaxytrue.d, omaxytrue.ull, omaxulperr,
           (double)omaxiulperr, 100.0 * omaxncorrnd / ntests);

    fprintf(fptr,"Min \n");

    __int64_t ominiulperr = ominytest.i - ominytrue.i;
    fprintf(fptr,"%3d %2d %5d %5s %10.2e %08lX %10.2e %016lX %10.2e %016lX %10.2e "
           "%10.2e %6.2f\n",
           minid, 0, 0, "     ",
           mantissa_to_double(ominx), (__int64_t)ominx, ominytest.d,
           ominytest.ull, ominytrue.d, ominytrue.ull, ominulperr,
           (double)ominiulperr, 100.0 * ominncorrnd / ntests);
    fprintf(fptr,"\n");

    fprintf(fptr,"\n");
    fprintf(fptr,"Histogram of exponent %d\n", j);
    fprintf(fptr,"Number of test %lu\n", ontests);
    fprintf(fptr,"===neg to pos===     ===count===  =====absolute===== \n");
    fprintf(fptr,"[-inf,-10) %6.2f%%  %10lu\n", 100.0 * n_ulp_n10down / ontests,
           n_ulp_n10down);
    fprintf(fptr,"  [-10,-5) %6.2f%%  %10lu\n", 100.0 * n_ulp_n10 / ontests,
           n_ulp_n10);
    fprintf(fptr,"   [-5,-2) %6.2f%%  %10lu\n", 100.0 * n_ulp_n5 / ontests,
           n_ulp_n5);
    fprintf(fptr,"   [-2,-1) %6.2f%%  %10lu\n", 100.0 * n_ulp_n2 / ontests,
           n_ulp_n2);
    fprintf(fptr," [-1,-0.5) %6.2f%%  %10lu\n", 100.0 * n_ulp_n1 / ontests,
           n_ulp_n1);
    fprintf(fptr,"  (-0.5,0] %6.2f%%  %10lu\n", 100.0 * n_ulp_n05 / ontests,
           n_ulp_n05);
    fprintf(fptr,"   (0,0.5] %6.2f%%  %10lu     [0,0.5] %6.2f%%\n",
           100.0 * n_ulp_05 / ontests, n_ulp_05,
           100.0 * (n_ulp_05 + n_ulp_n05) / ontests);
    fprintf(fptr,"   (0.5,1] %6.2f%%  %10lu     (0.5,1] %6.2f%%\n",
           100.0 * n_ulp_1 / ontests, n_ulp_1,
           100.0 * (n_ulp_1 + n_ulp_n1) / ontests);
    fprintf(fptr,"     (1,2] %6.2f%%  %10lu       (1,2] %6.2f%%\n",
           100.0 * n_ulp_2 / ontests, n_ulp_2,
           100.0 * (n_ulp_2 + n_ulp_n2) / ontests);
    fprintf(fptr,"     (2,5] %6.2f%%  %10lu       (2,5] %6.2f%%\n",
           100.0 * n_ulp_5 / ontests, n_ulp_5,
           100.0 * (n_ulp_5 + n_ulp_n5) / ontests);
    fprintf(fptr,"    (5,10] %6.2f%%  %10lu      (5,10] %6.2f%%\n",
           100.0 * n_ulp_10 / ontests, n_ulp_10,
           100.0 * (n_ulp_10 + n_ulp_n10) / ontests);
    fprintf(fptr,"  (10,inf] %6.2f%%  %10lu    (10,inf] %6.2f%%\n",
           100.0 * n_ulp_10up / ontests, n_ulp_10up,
           100.0 * (n_ulp_10up + n_ulp_n10down) / ontests);
    fprintf(fptr,"\n");
    // end j
  }

  if (exp > 0) {
    fprintf(fptr,"\n");
    fprintf(fptr,"Histogram of exponent 0 to exponent %d\n", exp);
    fprintf(fptr,"Number of test %lu\n", oontests);
    fprintf(fptr,"****neg to pos****     *****absolute*****\n");
    fprintf(fptr,"[-inf,-10) %6.2f%%\n", 100.0 * on_ulp_n10down / oontests);
    fprintf(fptr,"  [-10,-5) %6.2f%%\n", 100.0 * on_ulp_n10 / oontests);
    fprintf(fptr,"   [-5,-2) %6.2f%%\n", 100.0 * on_ulp_n5 / oontests);
    fprintf(fptr,"   [-2,-1) %6.2f%%\n", 100.0 * on_ulp_n2 / oontests);
    fprintf(fptr," [-1,-0.5) %6.2f%%\n", 100.0 * on_ulp_n1 / oontests);
    fprintf(fptr,"  (-0.5,0] %6.2f%%\n", 100.0 * on_ulp_n05 / oontests);
    fprintf(fptr,"   (0,0.5] %6.2f%%  [0,0.5] %6.2f%%\n",
           100.0 * on_ulp_05 / oontests,
           100.0 * (on_ulp_05 + on_ulp_n05) / oontests);
    fprintf(fptr,"   (0.5,1] %6.2f%%  (0.5,1] %6.2f%%\n", 100.0 * on_ulp_1 / oontests,
           100.0 * (on_ulp_1 + on_ulp_n1) / oontests);
    fprintf(fptr,"     (1,2] %6.2f%%    (1,2] %6.2f%%\n", 100.0 * on_ulp_2 / oontests,
           100.0 * (on_ulp_2 + on_ulp_n2) / oontests);
    fprintf(fptr,"     (2,5] %6.2f%%    (2,5] %6.2f%%\n", 100.0 * on_ulp_5 / oontests,
           100.0 * (on_ulp_5 + on_ulp_n5) / oontests);
    fprintf(fptr,"    (5,10] %6.2f%%   (5,10] %6.2f%%\n",
           100.0 * on_ulp_10 / oontests,
           100.0 * (on_ulp_10 + on_ulp_n10) / oontests);
    fprintf(fptr,"  (10,inf] %6.2f%% (10,inf] %6.2f%%\n",
           100.0 * on_ulp_10up / oontests,
           100.0 * (on_ulp_10up + on_ulp_n10down) / oontests);
    fprintf(fptr,"\n");
  }

  printf("\n===================================================\n");
  printf("RECIP SINGLE PRECISION LOG REPORT:\n");
  printf("Results for recip single precision recorded in:\n\tlog/recip_sp.log\n");
  printf("Ulp error interval:\n\t [%.5f, %.5f]\n", ominulperr, omaxulperr);
  printf("===================================================\n\n");
  sleep(2);
  
  // end
}

int main(int argc, const char *argv[]) {

  // fprintf(fptr,"Test rounding %s \n",rounding);

  __uint64_t T_Table[RECIP_TABLE_SIZE];
  __uint64_t S_Table[RECIP_TABLE_SIZE];
  __uint64_t *qTable;
  __uint8_t A_Table[RECIP_TABLE_SIZE];
  FILE *fptr;

  fptr = fopen("../../../log/recip_sp.log", "w");

  // Allocate heap space for qTable
  qTable = malloc(sizeof(__uint64_t) * QTableSize);

  makeTTable(T_Table);
  makeSTable(S_Table, T_Table);
  makeQTable(qTable, T_Table, S_Table);
  makeATable(A_Table, qTable);

  test(0, 0, 255 /*RECIP_TABLE_SIZE-1*/, T_Table, S_Table, A_Table, fptr);

  fprintf(fptr,"Testing sign 1 bits, exp 8 bits, significand %d bits \n", fracBits);
  fprintf(fptr,"Testing reciprocal table : %d bits table, size %d \n", fracBitsTable,
         RECIP_TABLE_SIZE);
  fprintf(fptr,"Testing input(highf) %d bits \n", highf);
  fprintf(fptr,"Testing lowu %d bits \n", lowu);

#ifdef SINGLE
  fprintf(fptr,"Testing single precision \n");
#ifdef RU
  fprintf(fptr,"Testing round to positive infinity \n");
#endif
#ifdef RD
  fprintf(fptr,"Testing round to negative infinity \n");
#endif
#ifdef RZ
  fprintf(fptr,"Testing round to zero \n");
#endif
#ifdef RN
  fprintf(fptr,"Testing round to nearest even \n");
#endif

#else
  fprintf(fptr,"Testing double precision \n");
#endif

  free(qTable);
  fclose(fptr);
  return 0;
}
