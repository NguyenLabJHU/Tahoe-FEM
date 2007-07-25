/* $Id: SuperLUMatrixT.cpp,v 1.11 2005-08-01 03:26:30 paklein Exp $ */
#include "SuperLUMatrixT.h"

/* library support */
#ifdef __SUPERLU__

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
SuperLUMatrixT::SuperLUMatrixT(ostream& out, int check_code, bool symmetric, bool print_stat, 
	IterRefine_t refine, const CommunicatorT& comm):
	GlobalMatrixT(out, check_code, comm),
	fIsSymFactorized(false),
	fIsNumFactorized(false)
{
	const char caller[] = "SuperLUMatrixT::SuperLUMatrixT";

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

    /* Set the default input options:
		options.Fact = DOFACT;
		options.Equil = YES;
    	options.ColPerm = COLAMD;
		options.DiagPivotThresh = 1.0;
    	options.Trans = NOTRANS;
    	options.IterRefine = NOREFINE;
    	options.SymmetricMode = NO;
    	options.PivotGrowth = NO;
    	options.ConditionNumber = NO;
    	options.PrintStat = YES; */
    set_default_options(&foptions);
    foptions.PrintStat = (print_stat) ? YES : NO;
    
#if __option (extended_errorcheck)
    foptions.PrintStat = YES;
#else
    foptions.PrintStat = NO;
#endif
	foptions.SymmetricMode = (symmetric) ? YES : NO;		
}

/* Destructor */	
SuperLUMatrixT::~SuperLUMatrixT(void)
{
	/* free the matrix */
	Destroy_CompCol_Matrix(&fA);
        Destroy_Dense_Matrix(&fX);
        free(fB.Store);
	/* free upper and lower factors */
	if (fIsSymFactorized) { 
		Destroy_SuperNode_Matrix(&fL);
		Destroy_CompCol_Matrix(&fU);
	} 
}

/* set the internal matrix structure.
* NOTE: do not call Initialize() until equation topology has been set
* with AddEquationSet() for all equation sets */
void SuperLUMatrixT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	const char caller[] = "SuperLUMatrixT::Initialize";

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
	fetree.Dimension(fLocNumEQ);

	/* structure could be changing, so get rid of old factors etc. */
	if (fIsSymFactorized) {
		Destroy_SuperNode_Matrix(&fL);
		fL.nrow = 0;
		fL.ncol = 0;
		fL.Store = NULL;
		Destroy_CompCol_Matrix(&fU);
		fU.nrow = 0;
		fU.ncol = 0;
		fU.Store = NULL;
		fIsSymFactorized = false;
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
	fequed = 'N';
}

/* write information to output stream after SuperLUMatrixT::Initialize
 * has been called */
void SuperLUMatrixT::Info(ostream& out)
{
	/* inherited */
	GlobalMatrixT::Info(out);

	/* output */
	NCformat *A = (NCformat*) fA.Store;	
	out <<" Number of nonzeros in global matrix = " << A->nnz <<"\n"<<endl;
}

/* set all matrix values to 0.0 */
void SuperLUMatrixT::Clear(void)
{
	/* inherited */
	GlobalMatrixT::Clear();

	/* clear entries */
	NCformat *A = (NCformat*) fA.Store;
	memset(A->nzval, 0, sizeof(double)*A->colptr[fLocNumEQ]);
	
	/* no equilibration */
	fequed = 'N';
	fIsNumFactorized = false;
}

/* add element group equations to the overall topology.
 * NOTE: assembly positions (equation numbers) = 1...fLocNumEQ
 * equations can be of fixed size (iArray2DT) or
 * variable length (RaggedArray2DT) */
void SuperLUMatrixT::AddEquationSet(const iArray2DT& eqset)
{
	fEqnos.AppendUnique (&eqset);
}

/* see AddEquationSet above */
void SuperLUMatrixT::AddEquationSet(const RaggedArray2DT<int>& eqset)
{
	fRaggedEqnos.AppendUnique (&eqset);
}

/* assemble the element contribution into the LHS matrix - assumes
* that elMat is square (n x n) and that eqnos is also length n.
* NOTE: assembly positions (equation numbers) = 1...fLocNumEQ */
void SuperLUMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos)
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

void SuperLUMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& row_eqnos,
	const ArrayT<int>& col_eqnos)
{
#pragma unused(elMat)
#pragma unused(row_eqnos)
#pragma unused(col_eqnos)
	ExceptionT::GeneralFail("SuperLUMatrixT::Assemble", "non-square not implemented");
}

void SuperLUMatrixT::Assemble(const nArrayT<double>& diagonal_elMat, const ArrayT<int>& eqnos)
{
#pragma unused(diagonal_elMat)
#pragma unused(eqnos)
	ExceptionT::GeneralFail("SuperLUMatrixT::Assemble", "diagonal not implemented");
}

/* number scope and reordering */
GlobalMatrixT::EquationNumberScopeT SuperLUMatrixT::EquationNumberScope(void) const
{
	return kLocal;
}

bool SuperLUMatrixT::RenumberEquations(void) const { return false; }

SuperLUMatrixT& SuperLUMatrixT::operator=(const SuperLUMatrixT&)
{
	ExceptionT::GeneralFail("SuperLUMatrixT::operator=", "not implemented");
	return *this;
}

/* return a clone of self */
GlobalMatrixT* SuperLUMatrixT::Clone(void) const {

	ExceptionT::GeneralFail("SuperLUMatrixT::operator=", "not implemented");
	return NULL;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* solution driver */
void SuperLUMatrixT::BackSubstitute(dArrayT& result)
{
	const char caller[] = "SuperLUMatrixT::BackSubstitute";
	
	/* needs symbolic and numeric factorization */
	if (!fIsSymFactorized)
		foptions.Fact = DOFACT;
	else if (!fIsNumFactorized) /* compute numeric factorization assuming same sparsity */
		foptions.Fact = SamePattern_SameRowPerm;
	else /* just solve linear system */
		foptions.Fact = FACTORED;

	/* rhs into B */
	fB.nrow = fLocNumEQ;
	DNformat* BStore = (DNformat*) fB.Store;
	BStore->lda   = fLocNumEQ;
	BStore->nzval = result.Pointer();

    /* Initialize the statistics variables. */
	SuperLUStat_t stat;    
    StatInit(&stat);
    
	/* call SuperLU */
	int info;
	mem_usage_t mem_usage;
	int lwork = 0; /* allocate space internally */
	void* work = NULL;
	double recip_pivot_growth;
	double rcond;
	double ferr;
	double berr;

	dgssvx(&foptions, &fA, fperm_c.Pointer(), fperm_r.Pointer(), fetree.Pointer(), &fequed,
		fR.Pointer(), fC.Pointer(), &fL, &fU, work, lwork,
		&fB, &fX, &recip_pivot_growth, &rcond, &ferr, &berr, &mem_usage, &stat, &info);

	/* check results */
	if (info != 0)
		ExceptionT::BadJacobianDet(caller, "dgssvx return %d with estimated condition number %g", info, rcond);

	/* report statistics */
    if (foptions.PrintStat) StatPrint(&stat);
    StatFree(&stat);

	/* always fully factorized on exit */
	fIsSymFactorized = true;
	fIsNumFactorized = true;

	/* copy result */
	DNformat* XStore = (DNformat*) fX.Store;
	result = (double*) XStore->nzval;
}

/* check functions */
void SuperLUMatrixT::PrintAllPivots(void) const
{
// temp - not implemented yet. Maybe inappropriate
}

void SuperLUMatrixT::PrintZeroPivots(void) const
{
// temp - not implemented yet. Maybe inappropriate
}

void SuperLUMatrixT::PrintLHS(bool force) const
{
	if (!force || fCheckCode != GlobalMatrixT::kPrintLHS)
		return;

	fOut << "\nLHS matrix:\n\n";
	fOut << (*this) << "\n\n";
}

/* (over)estimate the number of nonzeros, based on the equation sets */
void SuperLUMatrixT::EstimateNNZ (int *colLength, int &totalnnz)
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
void SuperLUMatrixT::InsertEquations (NCformat *A, int *colLength, int &nnz)
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
void SuperLUMatrixT::InsertNZ (NCformat *A, int *colLength, int &nnz, int c,
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
void SuperLUMatrixT::CompressColumns (NCformat *A, const int *colLength)
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

#endif /* __SUPERLU__ */
