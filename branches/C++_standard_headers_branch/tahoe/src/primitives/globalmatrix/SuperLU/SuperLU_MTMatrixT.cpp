/* $Id: SuperLU_MTMatrixT.cpp,v 1.3 2005-05-29 04:24:30 paklein Exp $ */
#include "SuperLU_MTMatrixT.h"

/* library support */
#ifdef __SUPERLU_MT__

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include "toolboxConstants.h"
#include "ExceptionT.h"

/* types that we use in these methods */
#include "dMatrixT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "RaggedArray2DT.h"
#include "ElementMatrixT.h"

using namespace Tahoe;

/* constructor */
SuperLU_MTMatrixT::SuperLU_MTMatrixT(ostream& out, int check_code, int num_threads, const CommunicatorT& comm):
	GlobalMatrixT(out, check_code, comm),
	fNumThreads(num_threads),
	fIsSymFactorized(false),
	fIsNumFactorized(false)
{
	const char caller[] = "SuperLU_MTMatrixT::SuperLU_MTMatrixT";

	/* set up NULL structures */
	fA.Stype = SLU_NC; /* column-wise, no supernodes */
	fA.Dtype = SLU_D;  /* storing doubles */
	fA.Mtype = SLU_GE; /* general matrix */
	fA.nrow = 0;
	fA.ncol = 0;
	fA.Store = malloc(sizeof(NCformat));
	if (!fA.Store) ExceptionT::OutOfMemory(caller);

	NCformat *A = (NCformat*) fA.Store;
	A->nnz = 0;
	A->nzval = NULL;
	A->rowind = NULL;
	A->colptr = NULL;

	/* The only important thing to initialize in L and U are pointers */
	fL.Store = NULL;
	fU.Store = NULL;

	/* rhs vector */
	fB.Stype = SLU_DN;
	fB.Dtype = SLU_D;
	fB.Mtype = SLU_GE;
	fB.nrow = 0;
	fB.ncol = 1;
	fB.Store = malloc(sizeof(DNformat));
	if (!fB.Store) ExceptionT::OutOfMemory(caller);

	DNformat* BStore = (DNformat*) fB.Store;
	BStore->lda = 0;
	BStore->nzval = NULL;

	/* solution  vector */
	fX.Stype = SLU_DN;
	fX.Dtype = SLU_D;
	fX.Mtype = SLU_GE;
	fX.nrow = 0;
	fX.ncol = 1;
	fX.Store = malloc(sizeof(DNformat));
	if (!fX.Store) ExceptionT::OutOfMemory(caller);

	DNformat* XStore = (DNformat*) fX.Store;
	XStore->lda = 0;
	XStore->nzval = NULL;

    /* set sefault parameters to control factorization. */
    foptions.nprocs       = fNumThreads;
    foptions.fact         = DOFACT;
    foptions.trans        = NOTRANS;
    foptions.refact       = NO;
    foptions.panel_size   = sp_ienv(1);
    foptions.relax        = sp_ienv(2);
    foptions.diag_pivot_thresh = 1.0;
    foptions.usepr        = NO;
    foptions.drop_tol     = 0.0;
    foptions.perm_c       = NULL;
    foptions.perm_r       = NULL;
    foptions.work         = NULL;
    foptions.lwork        = 0;
	foptions.etree        = NULL;
	foptions.colcnt_h     = NULL;
	foptions.part_super_h = NULL;

//#if __option (extended_errorcheck)
//    foptions.PrintStat = YES;
//#else
//    foptions.PrintStat = NO;
//#endif
//	foptions.SymmetricMode = (symmetric) ? YES : NO;		
}

/* Destructor */	
SuperLU_MTMatrixT::~SuperLU_MTMatrixT(void)
{
	/* free the matrix */
	Destroy_CompCol_Matrix(&fA);
	free(((DNformat*)fX.Store)->nzval);
    Destroy_SuperMatrix_Store(&fX);
    Destroy_SuperMatrix_Store(&fB);

	/* free upper and lower factors */
	if (fIsNumFactorized) {
		Destroy_SuperNode_SCP(&fL);
		Destroy_CompCol_NCP(&fU);
	}

	/* free workspace */
	if (foptions.etree) {
		free(foptions.etree);
		foptions.etree = NULL;
	}
	if (foptions.colcnt_h) {
		free(foptions.colcnt_h);
		foptions.colcnt_h = NULL;
	}
	if (foptions.part_super_h) {
		free(foptions.part_super_h);
		foptions.part_super_h = NULL;
	}
}

/* set the internal matrix structure.
* NOTE: do not call Initialize() until equation topology has been set
* with AddEquationSet() for all equation sets */
void SuperLU_MTMatrixT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	const char caller[] = "SuperLU_MTMatrixT::Initialize";

	/* inherited */
	GlobalMatrixT::Initialize(tot_num_eq, loc_num_eq, start_eq);

	/* check */
	if (tot_num_eq != loc_num_eq)
		ExceptionT::GeneralFail(caller, "tot eq %d != loc eq %d", tot_num_eq, loc_num_eq);

	/* A note on memory allocation: since SuperLU is a C library, */
	/* I use malloc/free instead of new/delete for the structures */
	/* that SuperLU accesses, just in case. */	

	/* solution vector */
	fX.nrow = fLocNumEQ;
	DNformat* XStore = (DNformat*) fX.Store;
	XStore->lda = fLocNumEQ;
	free(XStore->nzval);
	XStore->nzval = (double*) malloc(fLocNumEQ*sizeof(double));
	if (!XStore->nzval) ExceptionT::OutOfMemory(caller);

	/* dimension work space */
	fperm_c.Dimension(fLocNumEQ);
	fperm_r.Dimension(fLocNumEQ);

	/* structure could be changing, so get rid of old factors etc. */
	if (fIsNumFactorized) {
		Destroy_SuperNode_Matrix(&fL);
		fL.nrow = 0;
		fL.ncol = 0;
		fL.Store = NULL;
		Destroy_CompCol_Matrix(&fU);
		fU.nrow = 0;
		fU.ncol = 0;
		fU.Store = NULL;
		fIsNumFactorized = false;
	}

	/* configure A */
	fA.nrow = fLocNumEQ;
	fA.ncol = fLocNumEQ;
	NCformat *A = (NCformat*) fA.Store;

	free(A->colptr);
	A->colptr = (int*) calloc(fLocNumEQ+1, sizeof(int));
	if (!A->colptr) ExceptionT::OutOfMemory(caller);

	/* We now construct the sparsity pattern of A from the equation sets */
	/* check if A is already allocated */
	if (A->rowind) {
		free(A->rowind);
		free(A->nzval);
	}

	/* Begin by (over-)estimating the number of nonzeros per column */
	int *colLength = (int*) malloc (fLocNumEQ*sizeof(int));
	if (!colLength) ExceptionT::OutOfMemory(caller);
	EstimateNNZ(colLength, A->nnz);

	/* Now allocate enough room for row indices (wait until later for */
	/* the nonzero values themselves) */
	A->rowind = (int*) malloc (A->nnz*sizeof(int));
	if (!A->rowind) ExceptionT::OutOfMemory(caller);

	/* Using the upper bounds in colLength, set up provisional column */
	/* pointers */
	A->colptr[0] = 0;
	for (int i = 0; i < fLocNumEQ; ++i)
		A->colptr[i+1] = A->colptr[i] + colLength[i];

	/* and now go through all the elements, inserting all the equations */
	InsertEquations (A, colLength, A->nnz);

	/* Then we can compress A to eliminate spaces between columns */
	CompressColumns (A, colLength);
	free(colLength);  // no longer needed

	/* and finish by reallocating A->rowind and allocating A->nzval */
	A->rowind = (int*) realloc (A->rowind, A->nnz*sizeof(int));
	if (!A->rowind) ExceptionT::OutOfMemory(caller);

	A->nzval = (void*) malloc (A->nnz*sizeof(double));
	if (!A->nzval) ExceptionT::OutOfMemory(caller);

	/* scalings */
	fR.Dimension(fA.nrow);
	fC.Dimension(fA.ncol);

	/* clear stored equation sets */
	fEqnos.Clear();
	fRaggedEqnos.Clear();	

	/* reset flags/options */
	fIsSymFactorized = false;
	fequed = NOEQUIL;

	/* reset workspace */
	if (foptions.etree) {
		free(foptions.etree);
		foptions.etree = NULL;
	}
	if (foptions.colcnt_h) {
		free(foptions.colcnt_h);
		foptions.colcnt_h = NULL;
	}
	if (foptions.part_super_h) {
		free(foptions.part_super_h);
		foptions.part_super_h = NULL;
	}
}

/* write information to output stream after SuperLU_MTMatrixT::Initialize
 * has been called */
void SuperLU_MTMatrixT::Info(ostream& out)
{
	/* inherited */
	GlobalMatrixT::Info(out);

	/* output */
	NCformat *A = (NCformat*) fA.Store;	
	out <<" Number of nonzeros in global matrix = " << A->nnz <<"\n"<<endl;
}

/* set all matrix values to 0.0 */
void SuperLU_MTMatrixT::Clear(void)
{
	/* inherited */
	GlobalMatrixT::Clear();

	/* clear entries */
	NCformat *A = (NCformat*) fA.Store;
	memset(A->nzval, 0, sizeof(double)*A->colptr[fLocNumEQ]);
	
	/* no equilibration */
	fequed = NOEQUIL;
	fIsNumFactorized = false;
}

/* add element group equations to the overall topology.
 * NOTE: assembly positions (equation numbers) = 1...fLocNumEQ
 * equations can be of fixed size (iArray2DT) or
 * variable length (RaggedArray2DT) */
void SuperLU_MTMatrixT::AddEquationSet(const iArray2DT& eqset)
{
	fEqnos.AppendUnique (&eqset);
}

/* see AddEquationSet above */
void SuperLU_MTMatrixT::AddEquationSet(const RaggedArray2DT<int>& eqset)
{
	fRaggedEqnos.AppendUnique (&eqset);
}

/* assemble the element contribution into the LHS matrix - assumes
* that elMat is square (n x n) and that eqnos is also length n.
* NOTE: assembly positions (equation numbers) = 1...fLocNumEQ */
void SuperLU_MTMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos)
{
	/* element matrix format */
	ElementMatrixT::FormatT format = elMat.Format();

	/* two cases: element matrix is diagonal, or it's not. */
	if (format == ElementMatrixT::kDiagonal)
	{
		/* less work to do! We only add diagonal entries */
		const double *pelMat = elMat.Pointer();
		int inc = elMat.Rows() + 1; // how far apart diag entries are

		int nee = eqnos.Length();
		for (int eqdex = 0; eqdex < nee; ++eqdex)
		{
			int eqno = eqnos[eqdex] - 1;
			if (eqno > -1)   // active dof?
				(*this)(eqno,eqno) += *pelMat;
			pelMat += inc;
		}
	}
	else    /* otherwise there's a full matrix to deal with */
	{
		/* If it's symmetric and just a triangle is stored, */
		/* copy it over to get full storage */
		if (format == ElementMatrixT::kSymmetricUpper)
			elMat.CopySymmetric();

		int nee = eqnos.Length();  // number of equations for element
		for (int col = 0; col < nee; ++col)
		{
			int ceqno = eqnos[col] - 1;
			if (ceqno > -1)   // active dof?
			{
				for (int row = 0; row < nee; ++row)
				{
					int reqno = eqnos[row] - 1;
					if (reqno > -1) // active dof?
						(*this)(reqno,ceqno) += elMat(row,col);
				}
			}
		}
	}
}

void SuperLU_MTMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& row_eqnos,
	const ArrayT<int>& col_eqnos)
{
#pragma unused(elMat)
#pragma unused(row_eqnos)
#pragma unused(col_eqnos)
	ExceptionT::GeneralFail("SuperLU_MTMatrixT::Assemble", "non-square not implemented");
}

void SuperLU_MTMatrixT::Assemble(const nArrayT<double>& diagonal_elMat, const ArrayT<int>& eqnos)
{
#pragma unused(diagonal_elMat)
#pragma unused(eqnos)
	ExceptionT::GeneralFail("SuperLU_MTMatrixT::Assemble", "diagonal not implemented");
}

/* number scope and reordering */
GlobalMatrixT::EquationNumberScopeT SuperLU_MTMatrixT::EquationNumberScope(void) const
{
	return kLocal;
}

bool SuperLU_MTMatrixT::RenumberEquations(void) const { return false; }

SuperLU_MTMatrixT& SuperLU_MTMatrixT::operator=(const SuperLU_MTMatrixT&)
{
	ExceptionT::GeneralFail("SuperLU_MTMatrixT::operator=", "not implemented");
	return *this;
}

/* return a clone of self */
GlobalMatrixT* SuperLU_MTMatrixT::Clone(void) const {
	ExceptionT::GeneralFail("SuperLU_MTMatrixT::Clone", "not implemented");
	return (GlobalMatrixT*) this;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* solution driver */
void SuperLU_MTMatrixT::BackSubstitute(dArrayT& result)
{
	const char caller[] = "SuperLU_MTMatrixT::BackSubstitute";

	/* rhs into B */
	fB.nrow = fLocNumEQ;
	DNformat* BStore = (DNformat*) fB.Store;
	BStore->lda   = fLocNumEQ;
	BStore->nzval = result.Pointer();
	
	/* needs symbolic and numeric factorization */
	if (!fIsSymFactorized)
	{
		/* set flags */
		foptions.fact = DOFACT;
		foptions.refact = NO;
		
	    /* Get column permutation vector perm_c[], according to permc_spec:
	     *   permc_spec = 0: natural ordering 
	     *   permc_spec = 1: minimum degree ordering on structure of A'*A
	     *   permc_spec = 2: minimum degree ordering on structure of A'+A
	     *   permc_spec = 3: approximate minimum degree for unsymmetric matrices */    	
		if (foptions.fact == DOFACT) {
	    	int permc_spec = 1;
	    	get_perm_c(permc_spec, &fA, fperm_c.Pointer());
	    }
	    foptions.perm_c = fperm_c.Pointer();
	    foptions.perm_r = fperm_r.Pointer();	
	}
	else if (!fIsNumFactorized) /* compute numeric factorization assuming same sparsity */
	{
		/* set flags */
		foptions.fact = DOFACT;
		foptions.refact = YES;
	}
	else /* just solve linear system */
	{
		/* set flags */
		foptions.fact = FACTORED;
		foptions.refact = YES;
	}

	/* call SuperLU */
	int info;
	superlu_memusage_t mem_usage;
	double recip_pivot_growth;
	double rcond;
	double ferr;
	double berr;
	pdgssvx(fNumThreads, &foptions, &fA, 
		fperm_c.Pointer(), fperm_r.Pointer(), &fequed, fR.Pointer(), fC.Pointer(),
		&fL, &fU, &fB, &fX, 
		&recip_pivot_growth, &rcond, &ferr, &berr, &mem_usage, &info);

	/* check results */
	if (info != 0)
		ExceptionT::BadJacobianDet(caller, "pdgssvx return %d with estimated condition number %g", info, rcond);

	/* always fully factorized on exit */
	fIsSymFactorized = true;
	fIsNumFactorized = true;

	/* copy result */
	DNformat* XStore = (DNformat*) fX.Store;
	result = (double*) XStore->nzval;
}

/* check functions */
void SuperLU_MTMatrixT::PrintAllPivots(void) const
{
// temp - not implemented yet. Maybe inappropriate
}

void SuperLU_MTMatrixT::PrintZeroPivots(void) const
{
// temp - not implemented yet. Maybe inappropriate
}

void SuperLU_MTMatrixT::PrintLHS(bool force) const
{
	if (!force || fCheckCode != GlobalMatrixT::kPrintLHS)
		return;

	fOut << "\nLHS matrix:\n\n";
	fOut << (*this) << "\n\n";
}

/* (over)estimate the number of nonzeros, based on the equation sets */
void SuperLU_MTMatrixT::EstimateNNZ (int *colLength, int &totalnnz)
{
	/* The estimate is simple: forget about overlap of nodes between */
	/* elements, and just add up all the nonzeros of all the element */
	/* stiffness matrices. */
	
	totalnnz = 0;
	memset (colLength, 0, fLocNumEQ*sizeof(int));

	const iArray2DT* peq;
	fEqnos.Top();
	while (fEqnos.Next(peq))
	{
		int nel = peq->MajorDim();  /* number of elements */
		int nee = peq->MinorDim();  /* number of element equations */
		for (int j = 0; j < nel; ++j)
		{
			const int *eleqnos = (*peq)(j);
			for (int k = 0; k < nee; ++k)
			{
				int eq = eleqnos[k]-1;
				if (eq > -1)
				{
					colLength[eq] += nee;
					totalnnz += nee;
				}
			}
		}
	}

	const RaggedArray2DT<int>* prageq;
	fRaggedEqnos.Top();
	while (fRaggedEqnos.Next(prageq))
	{
		int nel = prageq->MajorDim();  /* number of elements */
		for (int j = 0; j < nel; ++j)
			{
			int nee = prageq->MinorDim(j); /* no. element equations */
			const int *eleqnos = (*prageq)(j);
			for (int k = 0; k < nee; ++k)
			{
				int eq = eleqnos[k]-1;
				if (eq > -1)
				{
					colLength[eq] += nee;
					totalnnz += nee;
				}
			}
		}
	}
}

/* insert all the element equations into A */
void SuperLU_MTMatrixT::InsertEquations (NCformat *A, int *colLength, int &nnz)
{
	/* Reset the lengths */
	memset (colLength, 0, fLocNumEQ*sizeof(int));
	nnz = 0;

	const iArray2DT* peq;
	fEqnos.Top();
	while (fEqnos.Next(peq))
	{
		int nel = peq->MajorDim();  /* number of elements */
		int nee = peq->MinorDim();  /* number of element equations */
		for (int j = 0; j < nel; ++j)
		{
			const int *eleqnos = (*peq)(j);
			for (int k = 0; k < nee; ++k)
			{
				int eq = eleqnos[k]-1;
				if (eq > -1)
					InsertNZ (A, colLength, nnz, eq, nee, eleqnos);
			}
		}
	}

	const RaggedArray2DT<int>* prageq;
	fRaggedEqnos.Top();
	while (fRaggedEqnos.Next(prageq))
	{
		int nel = prageq->MajorDim();  /* number of elements */
		for (int j = 0; j < nel; ++j)
			{
			int nee = prageq->MinorDim(j); /* no. element equations */
			const int *eleqnos = (*prageq)(j);
			for (int k = 0; k < nee; ++k)
			{
				int eq = eleqnos[k]-1;
				if (eq > -1)
					InsertNZ (A, colLength, nnz, eq, nee, eleqnos);
			}
		}
	}
}

/* Insert the list nzlist of nonzeros (with length nzlen) into column c of
* matrix A, using colLengths to keep track of column lengths, and keeping
* nnz up-to-date. The columns of A will have nonzeros in ascending order.
* The entries of nzlist are FORTRAN-indexed (starting at 1, so entries 0 or
* less are ignored) whereas A and c are C-indexed (starting at 0).
*/
void SuperLU_MTMatrixT::InsertNZ (NCformat *A, int *colLength, int &nnz, int c,
int nzlen, const int *nzlist)
{
	for (int i = 0; i < nzlen; ++i)
	{
		int newnz = nzlist[i]-1;

		if (newnz < 0)
			continue;   /* ignore inactive dof's */

		/* insert newnz into column c */
		int j = 0;
		while (j < colLength[c])
		{
			int nzj = A->rowind[A->colptr[c]+j];
			if (newnz < nzj)
				break;
			else if (newnz == nzj)
				goto AlreadyThere; // column c already has it
			++j;
		}
		/* insert newnz before j */
		++colLength[c];  // new nonzero
		++nnz;
		/* move up list by one */
		for (int k = colLength[c]-1; k > j; --k)
			A->rowind[A->colptr[c]+k] = A->rowind[A->colptr[c]+k-1];
		/* put newnz in */
		A->rowind[A->colptr[c]+j] = newnz;

		AlreadyThere: {} // nothing to do if nzlist[i]-1 is there
	}
}

/* compress columns in A */
void SuperLU_MTMatrixT::CompressColumns (NCformat *A, const int *colLength)
{
	int writeto = 0,  /* where we are writing to in A->rowind */
	    colstart;     /* where a column will start in A->rowind */

	/* look over the columns */
	for (int i = 0; i < fLocNumEQ; ++i)
	{
		colstart = writeto;    /* where column i will start */
		/* copy i's list of nonzeros down to writeto */
		for (int j = A->colptr[i]; j < A->colptr[i]+colLength[i]; ++j)
		{
			A->rowind[writeto] = A->rowind[j];
			++writeto;
		}
		/* set column i's pointer to the new location in A->rowind */
		A->colptr[i] = colstart;
	}
	/* set the last pointer to writeto, which now should == nnz */
	A->colptr[fLocNumEQ] = writeto;
}

#endif /* __SUPERLU_MT__ */
