/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: md_timer_mpi.c,v $
 *
 * $Author: paklein $
 *
 * $Date: 2003-02-28 01:19:20 $
 *
 * $Revision: 1.1 $
 *
 * $Name: not supported by cvs2svn $
 *====================================================================*/

/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/
#include "az_aztec_defs.h"
#ifdef AZTEC_MPI
#include <mpi.h>
#endif

/* MPI timer */

extern double second(void);
double second(void)

{

#ifdef AZTEC_MPI
  return (MPI_Wtime());
#else 
  return (0.0);
#endif

} /* second */
