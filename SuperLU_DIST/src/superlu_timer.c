/* 
 * Purpose
 * ======= 
 *	Returns the time in seconds used by the process.
 *
 * Note: the timer function call is machine dependent. Use conditional
 *       compilation to choose the appropriate function.
 *
 */

#include "superlu_defs.h"

#ifdef SUN 
/*
 * 	It uses the system call gethrtime(3C), which is accurate to 
 *	nanoseconds. 
*/
#include <sys/time.h>
 
double SuperLU_timer_() {
    return ( (double)gethrtime() / 1e9 );
}

#elif defined ( UNIX_TIMER )

#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include <sys/time.h>

#ifndef CLK_TCK
#define CLK_TCK 60
#endif

double SuperLU_timer_()
{
    struct tms use;
    double tmp;
    times(&use);
    tmp = use.tms_utime;
    tmp += use.tms_stime;
    return (double)(tmp) / (double) CLK_TCK;
}

#else

#include <mpi.h>

double SuperLU_timer_()
{
    return MPI_Wtime();
}

#endif

