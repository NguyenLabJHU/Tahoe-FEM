/* $Id: LU_MT_driver.h,v 1.2 2005-04-18 05:45:53 paklein Exp $ */
#ifndef _LU_MT_DRIVER__H_
#define _LU_MT_DRIVER__H_

#ifdef __cplusplus
extern "C" {
#endif

/* SPOOLES headers */
#include "SPOOLESMT.h"

/* init data structures needed for multiples solves */
extern int LU_MT_driver_init(int matrix_type,
	int symmetry_flag, int pivoting_flag, int seed, int num_eq,
	int num_thread, void** ppLU_dat);

/* compute factorization */
extern int LU_MT_driver_factorize(int msg_lvl, const char* message_file,
	int num_entries, int* r, int* c, double* v, void** ppLU_dat);

/* back substitute - solve rhs */
extern int LU_MT_driver_solve(int msg_lvl, const char* message_file,
	double* rhs2out, void** ppLU_dat);

/* free data structures needed for multiples solves */
extern int LU_MT_driver_free(void** ppLU_dat);

#ifdef __cplusplus
}
#endif

#endif  /* _LU_MT_DRIVER__H_ */
