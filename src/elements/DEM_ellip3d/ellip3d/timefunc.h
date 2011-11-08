#ifndef TIMEFUNC_H
#define TIMEFUNC_H

#include <sys/time.h>

struct timeval gettimediff(const struct timeval &time1, const struct timeval &time2);
long int microseconds(const struct timeval &time1, const struct timeval &time2);
double seconds(const struct timeval &time1, const struct timeval &time2);

#endif
