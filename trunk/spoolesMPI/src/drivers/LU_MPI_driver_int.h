/* $Id: LU_MPI_driver_int.h,v 1.1 2004-09-07 06:41:58 paklein Exp $ */
#ifndef _LU_MPI_DRIVER_INT_H_
#define _LU_MPI_DRIVER_INT_H_

/* 'external' header */
#include "LU_MPI_driver.h"

#ifdef __cplusplus
extern "C" {
#endif

/* structure needed for repeated solves */
typedef struct _LU_MPI_driver_data {

	/* flags and dimensions */
	int matrix_type;
	int symmetry_flag;
	int pivoting_flag;
	int rand_seed;
	int num_eq;
	int num_row;
	
	/* data structures */
	FrontMtx*      frontmtx;
	SolveMap*      solvemap;
	IV*            vtxmapIV;	
	
	/* allocated during factorization */	
	ETree*         frontETree;
	IV*            oldToNewIV;
	IV*            newToOldIV;
	IV*            ownersIV;

} LU_MPI_driver_data;

#ifdef __cplusplus
}
#endif

#endif  /* _LU_MPI_DRIVER_INT_H_ */
