/* $Id: pspases_int.h,v 1.4 2005-01-15 05:37:45 paklein Exp $ */

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
#define myMPI_Waitall   MPI_Waitall

#endif

#endif /* PSPASES_INT_H */
