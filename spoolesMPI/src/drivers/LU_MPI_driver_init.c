/* $Id: LU_MPI_driver_init.c,v 1.1 2004-09-07 06:41:58 paklein Exp $ */
#include "LU_MPI_driver_int.h"

/* init data structures needed for multiples solves */
int LU_MPI_driver_init(int matrix_type,
	int symmetry_flag, int pivoting_flag, int rand_seed, int num_eq, int num_row,
	void** ppLU_dat)
{
	/* resolve pointer */
	LU_MPI_driver_data* pLU_dat = (LU_MPI_driver_data*)(*ppLU_dat);
	if (pLU_dat) /* should be NULL */
		return 0;
		
	/* allocate structure */
	pLU_dat = (LU_MPI_driver_data*) malloc(sizeof(LU_MPI_driver_data));
	*ppLU_dat = (void*) pLU_dat;

	/* set flags and dimensions */
	pLU_dat->matrix_type = matrix_type;
	pLU_dat->symmetry_flag = symmetry_flag;
	pLU_dat->pivoting_flag = pivoting_flag;
	pLU_dat->rand_seed = rand_seed;
	pLU_dat->num_eq = num_eq;
	pLU_dat->num_row = num_row;
	
	/* data structures computed during factorization */
	pLU_dat->frontmtx = NULL;
	pLU_dat->solvemap = NULL;
	pLU_dat->vtxmapIV = NULL;
	pLU_dat->frontETree = NULL;
	pLU_dat->ownersIV = NULL;
	pLU_dat->newToOldIV = NULL;
	pLU_dat->oldToNewIV = NULL;

	return 1;
}
