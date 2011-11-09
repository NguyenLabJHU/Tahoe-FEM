#ifndef TIMEFUNC_H
#define TIMEFUNC_H

#include <sys/time.h>

struct timeval timediff(const struct timeval &time1, const struct timeval &time2);
long int       timediffmsec(const struct timeval &time1, const struct timeval &time2);
double         timediffsec(const struct timeval &time1, const struct timeval &time2);

#endif
