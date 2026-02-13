#include "CHelpers.h"
#include <stdlib.h>

void c_srand48(long seedval) { srand48(seedval); }
long c_lrand48(void) { return lrand48(); }
