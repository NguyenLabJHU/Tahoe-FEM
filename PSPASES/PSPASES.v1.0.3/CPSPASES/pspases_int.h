/* $Id: pspases_int.h,v 1.3 2005-01-15 00:22:29 paklein Exp $ */

#ifndef PSPASES_INT_H
#define PSPASES_INT_H

#include "pspases_f2c.h"

/* not debugging */
#ifdef NDEBUG

/* skip BLAS wrappers */
#define mydsyrk_ dsyrk_

/* skip MPI wrappers */
#define myMPI_Isend     MPI_Isend
#define myMPI_Send      MPI_Send
#define myMPI_Get_count MPI_Get_count

#endif

#endif /* PSPASES_INT_H */
