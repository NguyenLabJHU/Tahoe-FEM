/* $Id: pspases_int.h,v 1.2 2005-01-05 16:51:31 paklein Exp $ */

#ifndef PSPASES_INT_H
#define PSPASES_INT_H

#include "pspases_f2c.h"

/* not debugging */
#ifdef NDEBUG

/* skip BLAS wrappers */
#define mydsyrk_ dsyrk_

/* skip MPI wrappers */
#define myMPI_Isend     MPI_Isend
#define myMPI_Get_count MPI_Get_count

#endif

#endif /* PSPASES_INT_H */
