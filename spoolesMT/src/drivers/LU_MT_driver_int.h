/* $Id: LU_MT_driver_int.h,v 1.1 2005-04-05 16:04:19 paklein Exp $ */
#ifndef _LU_MT_DRIVER_INT_H_
#define _LU_MT_DRIVER_INT_H_

/* 'external' header */
#include "SPOOLESMT.h"

#ifdef __cplusplus
extern "C" {
#endif

/* structure needed for repeated solves */
typedef struct _LU_MT_driver_data {

	/* flags and dimensions */
	int matrix_type;
	int symmetry_flag;
	int pivoting_flag;
	int rand_seed;
	int num_eq;
	int num_row;
	int n_thread;
	
	/* data structures */
	FrontMtx*      frontmtx;
	SolveMap*      solvemap;
	SubMtxManager* mtxmanager;
	
	/* allocated during factorization */	
	ETree*         frontETree;
	IV*            oldToNewIV;
	IV*            newToOldIV;
	IV*            ownersIV;

} _LU_MT_driver_data;

#ifdef __cplusplus
}
#endif

#endif  /* _LU_MT_DRIVER_INT_H_ */
