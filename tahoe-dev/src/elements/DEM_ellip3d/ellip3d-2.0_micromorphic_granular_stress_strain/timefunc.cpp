#include "timefunc.h"

namespace dem {

struct timeval timediff(const struct timeval &time1, const struct timeval &time2) {
    struct timeval diff;

    diff.tv_sec = time2.tv_sec - time1.tv_sec; 
    if( ( diff.tv_usec = time2.tv_usec - time1.tv_usec) < 0) {
	diff.tv_usec += 1000000;
	--diff.tv_sec;
    }

    return(diff); 
}

long int timediffmsec(const struct timeval &time1, const struct timeval &time2) {
    struct timeval diff = timediff(time1, time2);
    return(diff.tv_sec * 1000000 + diff.tv_usec); 
}

REAL timediffsec(const struct timeval &time1, const struct timeval &time2) {
    return( (REAL) timediffmsec(time1,time2) / 1.0e+6);
}

}

/*

1. calender time
#include <time.h>
time_t time(time_t *t);

#include <sys/time.h>
int gettimeofday(struct timeval * tv,struct timezone *tz);
struct timeval {
  time_t      tv_sec;     // seconds
  suseconds_t tv_usec;    // microseconds
};

#include <time.h>
char *asctime(const struct tm *tm);
char *asctime_r(const struct tm *tm, char *buf);
char *ctime(const time_t *timep);
char *ctime_r(const time_t *timep, char *buf);
struct tm *gmtime(const time_t *timep);
struct tm *gmtime_r(const time_t *timep, struct tm *result);
struct tm *localtime(const time_t *timep);
struct tm *localtime_r(const time_t *timep, struct tm *result);
time_t mktime(struct tm *tm);

struct tm {
int tm_sec;         // seconds
int tm_min;         // minutes
int tm_hour;        // hours
int tm_mday;        // day of the month
int tm_mon;         // month
int tm_year;        // year
int tm_wday;        // day of the week
int tm_yday;        // day in the year
int tm_isdst;       // daylight saving time
};

2. process time
command: time

#include <sys/times.h>
clock_t times(struct tms *buf);

struct tms {
clock_t tms_utime;  // user time
clock_t tms_stime;  // system time
clock_t tms_cutime; // user time of children
clock_t tms_cstime; // system time of children
};

*/

