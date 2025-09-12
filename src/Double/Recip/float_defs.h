#define INF (1.0/0.0)
#define ABS(x) ((x)>=0 ? (x) : (-(x)))

typedef union {float f; int i;} float_union_t;
typedef union {double d; unsigned long long ull;} double_union_t;

// structure with IEEE float information, but decomposed for convenience
typedef struct {
    int sign;
    int exp; // unbiased exponent, ie the power of 2
    signed long long significand; // fracBits bits, not including assumed 1.
} float_internal;


#define RECIP_TABLE_SIZE 256 // 512 = 2 ^ 9// (256) // 256 = 2 ^ 8
#define SQRT_TABLE_SIZE 512
#define QTableSize 1048576 //8192 //16384 //8192 //4096 = (2 ^ 4) * 256
#define ATableSize 256 
#define oneOverSize 0.00390625 // increment between sample points 1/256 = 0.00390625
// 1/512 = 0.001953125
// 1/1024 = 0.0009765625
#define bit15 (1ull<<15) /* long int*/
#define bit14 (1ull<<14)
#define bit23 (1ull<<23)
#define bit22 (1ull<<22)
#define bit31 (1ull<<31)
#define bit27 (1ull<<27)
#define bit29 (1ull<<29)
#define bit28 (1ull<<28)
#define bit26 (1ull<<26)

// 1.0 = 10(fracBits-1)'s0 0.999 = 01(fracBits-1)'s0
#define fracBits FRACBITS
// highf and lowu only work for values between [1.0,2.0)
#define highf 8
#define lowu ((fracBits) - highf)
#define fstExp (1ull<<fracBits)
#define savedBits SB

#define qTableBit 10
#define aTableBit 8
#define fracBitsTable FRACTTABLE
#define lowuTable (fracBitsTable - highf)
// SEEEEEEEEffffffffuuuuuuuuuuuuuuuu
