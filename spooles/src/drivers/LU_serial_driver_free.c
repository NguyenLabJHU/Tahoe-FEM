/* $Id: LU_serial_driver_free.c,v 1.1 2004-09-07 06:40:19 paklein Exp $ */
#include "LU_serial_driver_int.h"

/* free data structures needed for multiples solves */
int LU_serial_driver_free(void** ppLU_dat)
{
	/* resolve pointer */
	LU_serial_driver_data* pLU_dat = (LU_serial_driver_data*)(*ppLU_dat);
	if (!pLU_dat) /* skip if NULL */
		return 0;

	/* set flags and dimensions */
	pLU_dat->seed = 0;
	pLU_dat->num_eq = 0;
	
	/* initialize data structures */
/*	InpMtx_free(pLU_dat->mtxA);
	pLU_dat->mtxA = NULL; */
/*	Graph_free(pLU_dat->graph);
	pLU_dat->graph = NULL; */
	FrontMtx_free(pLU_dat->frontmtx);
	pLU_dat->frontmtx = NULL;
	SubMtxManager_free(pLU_dat->mtxmanager);
	pLU_dat->mtxmanager = NULL;

	/* only non-NULL if matrix is factorized */
	if (pLU_dat->frontETree) {
		ETree_free(pLU_dat->frontETree);
		pLU_dat->frontETree = NULL;
	}
	if (pLU_dat->newToOldIV) {
		IV_free(pLU_dat->newToOldIV);
		pLU_dat->newToOldIV = NULL;
	}
	if (pLU_dat->oldToNewIV) {
		IV_free(pLU_dat->oldToNewIV);
		pLU_dat->oldToNewIV = NULL;
	}
/*	if (pLU_dat->symbfacIVL) {
		IVL_free(pLU_dat->symbfacIVL);
		pLU_dat->symbfacIVL = NULL;
	} */

	/* free structure */
	free(pLU_dat);
	*ppLU_dat = NULL;

	return 1;
}
