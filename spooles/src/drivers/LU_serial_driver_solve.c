/* $Id: LU_serial_driver_solve.c,v 1.1 2004-09-07 06:40:19 paklein Exp $ */
#include "LU_serial_driver_int.h"

/* back substitute - solve rhs */
int LU_serial_driver_solve(int msglvl, const char* message_file, 
	double* rhs2out, void** ppLU_dat)
{
/*
   ------------------------------------------------------------------
   all-in-one program to solve A X = Y

   (1) read in matrix entries for A and form InpMtx object
   (2) read in right hand side for Y entries and form DenseMtx object
   (3) form Graph object, order matrix and form front tree
   (4) get the permutation, permute the front tree, matrix 
       and right hand side and get the symbolic factorization
   (5) initialize the front matrix object to hold the factor matrices
   (6) compute the numeric factorization
   (7) post-process the factor matrices
   (8) compute the solution
   (9) permute the solution into the original ordering

   created -- 98jun04, cca
   ------------------------------------------------------------------
*/

/*--------------------------------------------------------------------*/
/* resolve pointer to the solver data */
LU_serial_driver_data* pLU_dat = (LU_serial_driver_data*)(*ppLU_dat);

/* re-used from solver data */
SubMtxManager* mtxmanager = NULL;
/* InpMtx*        mtxA = NULL; */
/* Graph*         graph = NULL; */
FrontMtx*      frontmtx = NULL;

/* ETree*         frontETree = NULL; */
IV*            newToOldIV = NULL;
IV*            oldToNewIV = NULL;
/* IVL*           symbfacIVL = NULL; */

/* flags and dimensions */
int num_eq = 0;
int matrix_type = 0;
int symmetry_flag = 0;
int seed = 0;
int pivoting_flag = 0;
/*--------------------------------------------------------------------*/

/* temporary work space */

DenseMtx        *mtxY;
DenseMtx        *mtxX;
/* Chv             *rootchv = NULL; */ /* returned by numeric factorization */
/* ChvManager      *chvmanager = NULL; */ /* used for numeric factorization */
double          droptol = 0.0, tau = 100. ;
double          cpus[10] ;
FILE            /* *inputFile,*/ *msgFile ;
int             error, ient, irow, jcol, jrhs, jrow, /* msglvl, ncol, */
                nedges, nent, neqns, nrhs /*, nrow,  pivotingflag,  seed, 
                symmetryflag, type */ ;
int             stats[20] ;
/* int             *newToOld, *oldToNew; */
/* IVL             *adjIVL; */
/*--------------------------------------------------------------------*/
/*
   --------------------
   open message file
   --------------------
*/
if ( (msgFile = fopen(message_file, "a")) == NULL )
{
   fprintf(stderr, "\n SPOOLES:LU_serial_driver: unable to open message file %s\n",
		message_file) ;
   return -1;
}
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
/*
   --------------------
   resolve data
   --------------------
*/

if (!pLU_dat) {
	fprintf(stderr, "\n SPOOLES:LU_serial_driver_solve: pLU_dat is NULL\n",
	message_file) ;
   	return -1;
}
else {

	/* re-used from solver data */
	mtxmanager = pLU_dat->mtxmanager;
/*	mtxA = pLU_dat->mtxA; */
/*	graph = pLU_dat->graph; */
	frontmtx = pLU_dat->frontmtx;

/* 	frontETree = pLU_dat->frontETree; */
	newToOldIV = pLU_dat->newToOldIV;
	oldToNewIV = pLU_dat->oldToNewIV;
/* 	symbfacIVL = pLU_dat->symbfacIVL; */

	/* flags and dimensions */
	num_eq = pLU_dat->num_eq;
	matrix_type = pLU_dat->matrix_type;
	symmetry_flag = pLU_dat->symmetry_flag;
	seed = pLU_dat->seed;
	pivoting_flag = pLU_dat->pivoting_flag;
}	

/*--------------------------------------------------------------------*/


/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   STEP 1: create the InpMtx object and load in
           {r, c, v} values for
   --------------------------------------------
*/
neqns = num_eq;
/* performed during factorization */
/*--------------------------------------------------------------------*/

/*
   -----------------------------------------
   STEP 2: read the right hand side matrix Y
   -----------------------------------------
*/
mtxY = DenseMtx_new();
nrhs = 1;
DenseMtx_init(mtxY, matrix_type, 0, 0, neqns, nrhs, 1, neqns);
DenseMtx_zero(mtxY) ;
if (matrix_type == SPOOLES_REAL)
{
	for (irow = 0 ; irow < neqns ; irow++) 
		DenseMtx_setRealEntry(mtxY, irow, 0, rhs2out[irow]);
		/* NOTE: there is no function to copy a vector into a column */
}
else
{
   fprintf(stderr, "\n SPOOLES:LU_serial_driver: real matrices only");
   return -1;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n rhs matrix in original ordering") ;
   DenseMtx_writeForHumanEye(mtxY, msgFile) ;
   fflush(msgFile) ;
}

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   STEP 3 : find a low-fill ordering
   (1) create the Graph object
   (2) order the graph using multiple minimum degree
   -------------------------------------------------
*/
/* performed during factorization */
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   STEP 4: get the permutation, permute the front tree,
           permute the matrix and right hand side, and 
           get the symbolic factorization
   ----------------------------------------------------
*/
/* operations on matrix performed during factorization */
DenseMtx_permuteRows(mtxY, oldToNewIV);
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n right hand side matrix after permutation") ;
   DenseMtx_writeForHumanEye(mtxY, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   STEP 5: initialize the front matrix object
   ------------------------------------------
*/
/* performed during factorization */
/*--------------------------------------------------------------------*/

/*
   -----------------------------------------
   STEP 6: compute the numeric factorization
   -----------------------------------------
*/
/* performed during factorization */
/*--------------------------------------------------------------------*/

/*
   --------------------------------------
   STEP 7: post-process the factorization
   --------------------------------------
*/
/* performed during factorization */
/*--------------------------------------------------------------------*/

/*
   -------------------------------
   STEP 8: solve the linear system
   -------------------------------
*/
mtxX = DenseMtx_new() ;
DenseMtx_init(mtxX, matrix_type, 0, 0, neqns, nrhs, 1, neqns) ;
DenseMtx_zero(mtxX) ;
FrontMtx_solve(frontmtx, mtxX, mtxY, mtxmanager, 
               cpus, msglvl, msgFile) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n solution matrix in new ordering") ;
   DenseMtx_writeForHumanEye(mtxX, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/

/*
   -------------------------------------------------------
   STEP 9: permute the solution into the original ordering
           and copy into return value
   -------------------------------------------------------
*/
DenseMtx_permuteRows(mtxX, newToOldIV);
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n solution matrix in original ordering") ;
   DenseMtx_writeForHumanEye(mtxX, msgFile) ;
   fflush(msgFile) ;
}
if (matrix_type == SPOOLES_REAL)
{
	for (irow = 0 ; irow < neqns ; irow++) 
		/* NOTE: there is no function to copy a column into a vector */
		if (DenseMtx_realEntry(mtxX, irow, 0, rhs2out + irow) != 1)
		{
			fprintf(msgFile, "\n error permuting solution to original ordering\n") ;
			return -1;
		}		
}
else
{
   fprintf(stderr, "\n SPOOLES:LU_serial_driver: real matrices only");
   return -1;
}

/*--------------------------------------------------------------------*/

/*
   -----------
   free memory
   -----------
*/
DenseMtx_free(mtxX) ;
DenseMtx_free(mtxY) ;

/*--------------------------------------------------------------------*/
fclose(msgFile);
return 1; 
}
/*--------------------------------------------------------------------*/
