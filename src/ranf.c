#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>

/* Fortran interface to random routine */

void ranset_(unsigned *seed) { static char state[256];  initstate(*seed, state, 256); }
double ranf_() { return(random()/2147483648.); }
