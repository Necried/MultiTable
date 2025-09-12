#include <stdlib.h>
void makeTTable(__uint64_t* table);
void makeSTable(__uint64_t *table, __uint64_t *reciptable);
void makeQTable(__uint64_t *table, __uint64_t *recipTable, __uint64_t *linearTable);
void makeATable (__uint8_t *table,__uint64_t* qTable);

double calculate_ytrue(unsigned __int128 x);
double calculate_ytest(unsigned __int128 x,
                       __uint64_t *recipTable, __uint64_t *linearTable,
                       __uint8_t *aTable);
int calculate_id(unsigned __int128 x);

unsigned __int128 three_table_procedure(unsigned __int128 x, __uint64_t *T_Table,
                             __uint64_t *S_Table, __uint8_t *A_Table);