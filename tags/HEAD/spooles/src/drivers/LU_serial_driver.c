/* SPOOLES headers */
#include "misc.h"
#include "FrontMtx.h"
#include "SymbFac.h"

/* serial driver provided by in SPOOLES documentation */
int LU_serial_driver(int msglvl, const char* message_file, int matrix_type,
	int symmetry_flag, int pivoting_flag, int seed, int num_eq, double* rhs2out, 
	int num_entries, int* r, int* c, double* v)
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
/*char            *matrixFileName, *rhsFileName ; */
DenseMtx        *mtxY, *mtxX ;
Chv             *rootchv ;
ChvManager      *chvmanager  ;
SubMtxManager   *mtxmanager  ;
FrontMtx        *frontmtx ;
InpMtx          *mtxA ;
double          droptol = 0.0, tau = 100. ;
double          cpus[10] ;
ETree           *frontETree ;
FILE            /* *inputFile,*/ *msgFile ;
Graph           *graph ;
int             error, ient, irow, jcol, jrhs, jrow, /* msglvl, ncol, */
                nedges, nent, neqns, nrhs /*, nrow,  pivotingflag,  seed, 
                symmetryflag, type */ ;
int             *newToOld, *oldToNew ;
int             stats[20] ;
IV              *newToOldIV, *oldToNewIV ;
IVL             *adjIVL, *symbfacIVL ;
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
/*
   --------------------------------------------
   STEP 1: create the InpMtx object and load in
           {r, c, v} values for
   --------------------------------------------
*/
neqns = num_eq;
nent  = num_entries;
mtxA = InpMtx_new() ;
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
graph = Graph_new() ;
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
DenseMtx_permuteRows(mtxY, oldToNewIV) ;
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
   fprintf(msgFile, "\n\n right hand side matrix after permutation") ;
   DenseMtx_writeForHumanEye(mtxY, msgFile) ;
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
frontmtx   = FrontMtx_new() ;
mtxmanager = SubMtxManager_new() ;
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
ChvManager_free(chvmanager) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n factor matrix") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   fflush(msgFile) ;
}
if ( rootchv != NULL ) {
   fprintf(msgFile, "\n\n matrix found to be singular\n") ;
   exit(-1) ;
}
if ( error >= 0 ) {
   fprintf(msgFile, "\n\n error encountered at front %d", error) ;
   exit(-1) ;
}
/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   STEP 7: post-process the factorization
   --------------------------------------
*/
FrontMtx_postProcess(frontmtx, msglvl, msgFile) ;
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
		DenseMtx_realEntry(mtxX, irow, 0, rhs2out + irow);
		/* NOTE: there is no function to copy a column into a vector */
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
FrontMtx_free(frontmtx) ;
DenseMtx_free(mtxX) ;
DenseMtx_free(mtxY) ;
IV_free(newToOldIV) ;
IV_free(oldToNewIV) ;
InpMtx_free(mtxA) ;
ETree_free(frontETree) ;
IVL_free(symbfacIVL) ;
SubMtxManager_free(mtxmanager) ;
Graph_free(graph) ;
/*--------------------------------------------------------------------*/
fclose(msgFile);
return 1; }
/*--------------------------------------------------------------------*/
