/* $Id: LU_serial_driver_int.h,v 1.1 2004-09-07 06:40:19 paklein Exp $ */
#ifndef _LU_SERIAL_DRIVER_INT_H_
#define _LU_SERIAL_DRIVER_INT_H_

#ifdef __cplusplus
extern "C" {
#endif

/* 'external' header */
#include "LU_serial_driver.h"

/* SPOOLES headers */
#include "misc.h"
#include "FrontMtx.h"
#include "SymbFac.h"

/* structure needed for repeated solves */
typedef struct _LU_serial_driver_data {

	/* flags and dimensions */
	int matrix_type;
	int symmetry_flag;
	int pivoting_flag;
	int seed;
	int num_eq;
	
	/* data structures */
	/* InpMtx*        mtxA; */
	/* Graph*         graph; */
	FrontMtx*      frontmtx;
	SubMtxManager* mtxmanager;

	/* allocated during factorize */	
	ETree*         frontETree;
	IV*            newToOldIV;
	IV*            oldToNewIV;
	/* IVL*           symbfacIVL; */

} LU_serial_driver_data;

#ifdef __cplusplus
}
#endif

#endif  /* _LU_SERIAL_DRIVER_INT_H_ */
