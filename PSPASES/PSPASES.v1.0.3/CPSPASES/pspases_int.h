/* $Id: pspases_int.h,v 1.5 2005-01-15 08:18:28 paklein Exp $ */

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

#else /* use debugging wrappers */

/* skip BLAS wrappers */
#define mydsyrk_ dsyrk_d

/* skip MPI wrappers */
#define myMPI_Isend     MPI_Isend_d
#define myMPI_Send      MPI_Send_d
#define myMPI_Get_count MPI_Get_count_d
#define myMPI_Waitall   MPI_Waitall_d

#endif

#endif /* PSPASES_INT_H */
