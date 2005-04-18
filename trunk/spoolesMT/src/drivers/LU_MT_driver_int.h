/* $Id: LU_MT_driver_int.h,v 1.2 2005-04-18 05:45:54 paklein Exp $ */
#ifndef _LU_MT_DRIVER_INT_H_
#define _LU_MT_DRIVER_INT_H_

#ifdef __cplusplus
extern "C" {
#endif

/* 'external' header */
#include "SPOOLESMT.h"

/* SPOOLES headers */
#include "MT.spoolesMT.h"
#include "misc.h"
#include "SymbFac.h"

/* structure needed for repeated solves */
typedef struct _LU_MT_driver_data {

	/* flags and dimensions */
	int matrix_type;
	int symmetry_flag;
	int pivoting_flag;
	int rand_seed;
	int num_eq;
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

} LU_MT_driver_data;

#ifdef __cplusplus
}
#endif

#endif  /* _LU_MT_DRIVER_INT_H_ */
