/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: md_timer_generic.c,v $
 *
 * $Author: paklein $
 *
 * $Date: 2001-01-30 20:59:13 $
 *
 * $Revision: 1.1.1.1 $
 *
 * $Name: not supported by cvs2svn $
 *====================================================================*/
#ifndef lint
static char *cvs_timergen_id =
  "$Id: md_timer_generic.c,v 1.1.1.1 2001-01-30 20:59:13 paklein Exp $";
#endif

/*
     Clock rolls negative when shifts to 32nd bit
     wrap time is about  = 2147 seconds = about 36 minutes
     clockp = double(1<<30)*4.*1.e-6 is the proper increment
     PARAMETER (CLOCKP = 4294.967296)
*/

#include <time.h>

double second();

double second()
/*
 * Returns system cpu and wall clock time in seconds
 */
{
  time_t  itime;
  clock_t num_ticks, clock();
  double  cpu, wall;

  num_ticks=clock();
  if(num_ticks==-1) cpu=-1.;
  else cpu = ((double) num_ticks)/((double) CLOCKS_PER_SEC);
  /* wall= ((double) time(&itime)); */

  return(cpu);
}
