#ifndef CHELPERS_H
#define CHELPERS_H

#include <zlib.h>

// Wrappers for POSIX functions unavailable in Swift
void c_srand48(long seedval);
long c_lrand48(void);

#endif
