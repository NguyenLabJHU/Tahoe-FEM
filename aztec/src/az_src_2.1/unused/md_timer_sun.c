/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: md_timer_sun.c,v $
 *
 * $Author: paklein $
 *
 * $Date: 2003-02-28 01:19:20 $
 *
 * $Revision: 1.1 $
 *
 * $Name: not supported by cvs2svn $
 *====================================================================*/
#ifndef lint
static char *cvs_timersun_id =
  "$Id: md_timer_sun.c,v 1.1 2003-02-28 01:19:20 paklein Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


#include <sys/time.h>
#include <sys/resource.h>

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
extern double second(void);
/* extern int getrusage(void);*/

double second(void)

{
  double mytime;                  /* elapsed time in seconds */

  struct rusage rusage;


  getrusage(RUSAGE_SELF, &rusage);

  mytime = ((rusage.ru_utime.tv_sec + rusage.ru_stime.tv_sec) +
          1.0e-6 * (rusage.ru_utime.tv_usec + rusage.ru_stime.tv_usec));

  return mytime;

} /* second */
