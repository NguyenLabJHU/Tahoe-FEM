#ifndef TIMEFUNC_H
#define TIMEFUNC_H

#include "realtypes.h"
#include <sys/time.h>

namespace dem {

struct timeval timedifsf(const struct timeval &time1, const struct timeval &time2);
long int       timediffmsec(const struct timeval &time1, const struct timeval &time2);
REAL           timediffsec(const struct timeval &time1, const struct timeval &time2);

}

#endif
