/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: md_timer_ncube.c,v $
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
static char *cvs_timerncube_id =
  "$Id: md_timer_ncube.c,v 1.1 2003-02-28 01:19:20 paklein Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


/* nCUBE timer */

extern double second(void);
extern long int n_time(int *);
double second(void)

{

  long int nstart;
  int      nstart_micro;

  nstart = n_time(&nstart_micro);

  return ((double) nstart + ((double) nstart_micro) * 1.0e-6);

} /* second */
