#ifndef _DOUBLE_RECIP_UTILS
#define _DOUBLE_RECIP_UTILS

#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x) do {} while (0)
#endif

#define STR_(x)  #x
#define STR(x) STR_(x) // Converts macro to string

#endif // _DOUBLE_RECIP_UTILS
