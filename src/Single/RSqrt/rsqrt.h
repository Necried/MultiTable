#include <stdlib.h>
void makeTTable(__uint64_t* table);
void makeSTable(__uint64_t *table, __uint64_t *T_Table);
void makeQTable(__uint64_t *table, __uint64_t *T_Table, __uint64_t *S_Table);
void makeATable (__uint8_t *table,__uint64_t* Q_Table);

double calculate_ytrue(unsigned __int128 x);
double calculate_ytest(unsigned __int128 x,
                       __uint64_t *T_Table, __uint64_t *S_Table,
                       __uint8_t *A_Table);
unsigned __int128 three_table_procedure(unsigned __int128 x,
                                        __uint64_t *T_Table,
                                        __uint64_t *S_Table,
                                        __uint8_t *A_Table);                       
int calculate_id(unsigned __int128 x);