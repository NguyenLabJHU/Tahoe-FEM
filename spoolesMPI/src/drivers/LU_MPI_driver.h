/* $Id: LU_MPI_driver.h,v 1.1 2004-09-07 06:41:58 paklein Exp $ */
#ifndef _LU_MPI_DRIVER__H_
#define _LU_MPI_DRIVER__H_

/* SPOOLES headers */
#include "MPI.spoolesMPI.h"
#include "timings.h"

#ifdef __cplusplus
extern "C" {
#endif

/* init data structures needed for multiples solves */
extern int LU_MPI_driver_init(int matrix_type,
	int symmetry_flag, int pivoting_flag, int seed, int num_eq, int num_row,
	void** ppLU_dat);

/* compute factorization */
extern int LU_MPI_driver_factorize(int msg_lvl, const char* message_file,
	int num_entries, int* r, int* c, double* v, void** ppLU_dat,
	MPI_Comm* comm);

/* back substitute - solve rhs */
extern int LU_MPI_driver_solve(int msg_lvl, const char* message_file,
	int* r_rhs, double* rhs2out, void** ppLU_dat,
	MPI_Comm* comm);

/* free data structures needed for multiples solves */
extern int LU_MPI_driver_free(void** ppLU_dat);

#ifdef __cplusplus
}
#endif

#endif  /* _LU_MPI_DRIVER__H_ */
