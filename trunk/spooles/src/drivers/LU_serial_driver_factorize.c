/* $Id: LU_serial_driver_factorize.c,v 1.1 2004-09-07 06:40:19 paklein Exp $ */
#include "LU_serial_driver_int.h"

/* compute factorization */
int LU_serial_driver_factorize(int msglvl, const char* message_file, 
	int num_entries, int* r, int* c, double* v,
	void** ppLU_dat)
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

ETree*         frontETree = NULL;
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

InpMtx*         mtxA = NULL;
Graph*          graph = NULL;
IVL*            symbfacIVL = NULL;
Chv             *rootchv = NULL; /* returned by numeric factorization */
ChvManager      *chvmanager = NULL; /* used for numeric factorization */
double          droptol = 0.0, tau = 100. ;
double          cpus[10] ;
FILE            /* *inputFile,*/ *msgFile ;
int             error, ient, irow, jcol, jrhs, jrow, /* msglvl, ncol, */
                nedges, nent, neqns, nrhs /*, nrow,  pivotingflag,  seed, 
                symmetryflag, type */ ;
int             stats[20] ;
int             *newToOld, *oldToNew;
IVL             *adjIVL;
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
	fprintf(stderr, "\n SPOOLES:LU_serial_driver_factorize: pLU_dat is NULL\n",
	message_file) ;
   	return -1;
}
else {

	/* structures allocated during init */
	mtxmanager = pLU_dat->mtxmanager;
/*	mtxA = pLU_dat->mtxA; */
/*	graph = pLU_dat->graph; */
	frontmtx = pLU_dat->frontmtx;

	/* flags and dimensions */
	num_eq = pLU_dat->num_eq;
	matrix_type = pLU_dat->matrix_type;
	symmetry_flag = pLU_dat->symmetry_flag;
	seed = pLU_dat->seed;
	pivoting_flag = pLU_dat->pivoting_flag;
}	

/*--------------------------------------------------------------------*/

/*
   --------------------------------------------
   STEP 1: create the InpMtx object and load in
           {r, c, v} values for
   --------------------------------------------
*/
neqns = num_eq;
nent  = num_entries;
mtxA = InpMtx_new();
InpMtx_init(mtxA, INPMTX_BY_ROWS, matrix_type, nent, neqns) ;
if (matrix_type == SPOOLES_REAL)
	/* enter all {row, column, value} triples */
	InpMtx_inputRealTriples(mtxA, num_entries, r, c, v);
else
{
   fprintf(stderr, "\n SPOOLES:LU_serial_driver: real matrices only");
   return -1;
}
InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n input matrix") ;
   InpMtx_writeForHumanEye(mtxA, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/

/* move to back substitute */
#if 0
/*
   -----------------------------------------
   STEP 2: read the right hand side matrix Y
   -----------------------------------------
*/
mtxY = DenseMtx_new(); in back substitute
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
#endif
/* move to back substitute */

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   STEP 3 : find a low-fill ordering
   (1) create the Graph object
   (2) order the graph using multiple minimum degree
   -------------------------------------------------
*/
graph = Graph_new();
adjIVL = InpMtx_fullAdjacency(mtxA) ;
nedges = IVL_tsize(adjIVL) ;
Graph_init2(graph, 0, neqns, 0, nedges, neqns, nedges, adjIVL,
            NULL, NULL) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n graph of the input matrix") ;
   Graph_writeForHumanEye(graph, msgFile) ;
   fflush(msgFile) ;
}
frontETree = orderViaMMD(graph, seed, msglvl, msgFile) ;
Graph_free(graph);
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front tree from ordering") ;
   ETree_writeForHumanEye(frontETree, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   STEP 4: get the permutation, permute the front tree,
           permute the matrix and right hand side, and 
           get the symbolic factorization
   ----------------------------------------------------
*/
oldToNewIV = ETree_oldToNewVtxPerm(frontETree) ;
oldToNew   = IV_entries(oldToNewIV) ;
newToOldIV = ETree_newToOldVtxPerm(frontETree) ;
newToOld   = IV_entries(newToOldIV) ;
ETree_permuteVertices(frontETree, oldToNewIV) ;
InpMtx_permute(mtxA, oldToNew, oldToNew) ;
if (  symmetry_flag == SPOOLES_SYMMETRIC
   || symmetry_flag == SPOOLES_HERMITIAN ) {
   InpMtx_mapToUpperTriangle(mtxA) ;
}
InpMtx_changeCoordType(mtxA, INPMTX_BY_CHEVRONS) ;
InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
/* DenseMtx_permuteRows(mtxY, oldToNewIV); move to back substitute */
symbfacIVL = SymbFac_initFromInpMtx(frontETree, mtxA) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n old-to-new permutation vector") ;
   IV_writeForHumanEye(oldToNewIV, msgFile) ;
   fprintf(msgFile, "\n\n new-to-old permutation vector") ;
   IV_writeForHumanEye(newToOldIV, msgFile) ;
   fprintf(msgFile, "\n\n front tree after permutation") ;
   ETree_writeForHumanEye(frontETree, msgFile) ;
   fprintf(msgFile, "\n\n input matrix after permutation") ;
   InpMtx_writeForHumanEye(mtxA, msgFile) ;
/*
   fprintf(msgFile, "\n\n right hand side matrix after permutation") ;
   DenseMtx_writeForHumanEye(mtxY, msgFile) ;
*/
   fprintf(msgFile, "\n\n symbolic factorization") ;
   IVL_writeForHumanEye(symbfacIVL, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   STEP 5: initialize the front matrix object
   ------------------------------------------
*/
/* frontmtx   = FrontMtx_new() ; re-used */
/* mtxmanager = SubMtxManager_new() ; re-used */
SubMtxManager_init(mtxmanager, NO_LOCK, 0) ;
FrontMtx_init(frontmtx, frontETree, symbfacIVL, matrix_type, symmetry_flag, 
              FRONTMTX_DENSE_FRONTS, pivoting_flag, NO_LOCK, 0, NULL, 
              mtxmanager, msglvl, msgFile) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   STEP 6: compute the numeric factorization
   -----------------------------------------
*/
chvmanager = ChvManager_new() ;
ChvManager_init(chvmanager, NO_LOCK, 1) ;
DVfill(10, cpus, 0.0) ;
IVfill(20, stats, 0) ;
rootchv = FrontMtx_factorInpMtx(frontmtx, mtxA, tau, droptol, 
             chvmanager, &error, cpus, stats, msglvl, msgFile) ;
ChvManager_free(chvmanager);
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n factor matrix") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   fflush(msgFile) ;
}
if ( rootchv != NULL ) {
   fprintf(msgFile, "\n\n matrix found to be singular\n") ;
	return -1;
}
if ( error >= 0 ) {
   fprintf(msgFile, "\n\n error encountered at front %d", error) ;
	return -1;
}
/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   STEP 7: post-process the factorization
   --------------------------------------
*/
FrontMtx_postProcess(frontmtx, msglvl, msgFile);
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n factor matrix after post-processing") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   -------------------------------
   STEP 8: solve the linear system
   -------------------------------
*/
/*
   -------------------------------------------------------
   STEP 9: permute the solution into the original ordering
           and copy into return value
   -------------------------------------------------------
*/

/* store reused data */
pLU_dat->frontETree = frontETree;
pLU_dat->newToOldIV = newToOldIV;
pLU_dat->oldToNewIV = oldToNewIV;
/* pLU_dat->symbfacIVL = symbfacIVL; */

/* clean up */
InpMtx_free(mtxA);
IVL_free(symbfacIVL); /* clears adjIVL as well */

fclose(msgFile);
return 1; 
}
/*--------------------------------------------------------------------*/
