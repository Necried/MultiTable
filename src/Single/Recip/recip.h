#include <stdlib.h>
void makeTTable(__uint64_t* table);
void makeSTable(__uint64_t *table, __uint64_t *reciptable);
void makeQTable(__uint64_t *table, __uint64_t *recipTable, __uint64_t *linearTable);
void makeATable (unsigned* table,__uint64_t* qTable);

double calculate_ytrue(unsigned __int128 x);
double calculate_ytest(unsigned __int128 x,
                       __uint64_t *recipTable, __uint64_t *linearTable,
                       unsigned *aTable);
int calculate_id(unsigned __int128 x);