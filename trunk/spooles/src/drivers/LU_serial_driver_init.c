/* $Id: LU_serial_driver_init.c,v 1.1 2004-09-07 06:40:19 paklein Exp $ */
#include "LU_serial_driver_int.h"

/* init data structures needed for multiples solves */
int LU_serial_driver_init(int matrix_type,
	int symmetry_flag, int pivoting_flag, int seed, int num_eq,
	void** ppLU_dat)
{
	/* resolve pointer */
	LU_serial_driver_data* pLU_dat = (LU_serial_driver_data*)(*ppLU_dat);
	if (pLU_dat) /* should be NULL */
		return 0;
		
	/* allocate structure */
	pLU_dat = (LU_serial_driver_data*) malloc(sizeof(LU_serial_driver_data));
	*ppLU_dat = (void*) pLU_dat;

	/* set flags and dimensions */
	pLU_dat->matrix_type = matrix_type;
	pLU_dat->symmetry_flag = symmetry_flag;
	pLU_dat->pivoting_flag = pivoting_flag;
	pLU_dat->seed = seed;
	pLU_dat->num_eq = num_eq;
	
	/* allocate data structures */
	pLU_dat->mtxmanager = SubMtxManager_new();
/*	pLU_dat->mtxA = InpMtx_new(); */
/*	pLU_dat->graph = Graph_new(); */
	pLU_dat->frontmtx = FrontMtx_new();

	/* allocated during factorization */
	pLU_dat->frontETree = NULL;
	pLU_dat->newToOldIV = NULL;
	pLU_dat->oldToNewIV = NULL;
/*	pLU_dat->symbfacIVL = NULL; */

	return 1;
}
