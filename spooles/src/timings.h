#ifndef _TIMINGS_
#define _TIMINGS_

#ifdef __MWERKS__
#include <time.h>
static clock_t TV;
#define MARKTIME(t) \
   TV = clock() ; \
   t = ((double) TV)/CLOCKS_PER_SEC;
#else
#include <sys/time.h>
static struct timeval  TV ;
static struct timezone TZ ;
#define MARKTIME(t) \
   gettimeofday(&TV, &TZ) ; \
   t = (TV.tv_sec + 0.000001*TV.tv_usec)
#endif
#endif /* _TIMINGS_ */
