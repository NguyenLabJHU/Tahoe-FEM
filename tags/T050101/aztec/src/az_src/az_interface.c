/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_interface.c,v $
 *
 * $Author: paklein $
 *
 * $Date: 2001-01-30 20:59:12 $
 *
 * $Revision: 1.1.1.1 $
 *
 * $Name: not supported by cvs2svn $
 *====================================================================*/
#ifndef lint
static char rcsid[] = "$Id: az_interface.c,v 1.1.1.1 2001-01-30 20:59:12 paklein Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/

#include "az_aztec.h"
#include "az_aztec_defs.h"

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_second(void)

{
  extern double second(void);

  return second();

} /* AZ_second */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_processor_info(int proc_config[])

{
  extern void get_parallel_info(int *proc, int *nprocs, int *dim);

  get_parallel_info(&(proc_config[AZ_node]), &(proc_config[AZ_N_procs]),
                    &(proc_config[AZ_dim]));

} /* AZ_processor_info */
