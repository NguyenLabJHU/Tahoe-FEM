/* $Id: LU_serial_driver.h,v 1.1 2004-09-07 06:40:19 paklein Exp $ */
#ifndef _LU_SERIAL_DRIVER__H_
#define _LU_SERIAL_DRIVER__H_

#ifdef __cplusplus
extern "C" {
#endif

/* SPOOLES definitions */
#include "SPOOLES.h"

/* init data structures needed for multiples solves */
extern int LU_serial_driver_init(int matrix_type,
	int symmetry_flag, int pivoting_flag, int seed, int num_eq,
	void** ppLU_dat);

/* compute factorization */
extern int LU_serial_driver_factorize(int msglvl, const char* message_file,
	int num_entries, int* r, int* c, double* v,
	void** ppLU_dat);

/* back substitute - solve rhs */
extern int LU_serial_driver_solve(int msglvl, const char* message_file, 
	double* rhs2out, void** ppLU_dat);

/* free data structures needed for multiples solves */
extern int LU_serial_driver_free(void** ppLU_dat);

#ifdef __cplusplus
}
#endif

#endif  /* _LU_SERIAL_DRIVER__H_ */
