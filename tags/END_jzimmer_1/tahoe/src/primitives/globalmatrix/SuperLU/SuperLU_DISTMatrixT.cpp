/* $Id: SuperLU_DISTMatrixT.cpp,v 1.1 2004-03-16 10:03:21 paklein Exp $ */
#include "SuperLU_DISTMatrixT.h"

/* library support options */
#ifdef __SUPERLU_DIST__
#ifdef __TAHOE_MPI__

#include "MSRBuilderT.h"
#include "dMatrixT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "ElementMatrixT.h"

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

using namespace Tahoe;

/* constructor */
SuperLU_DISTMatrixT::SuperLU_DISTMatrixT(ostream& out, int check_code, bool symmetric, CommunicatorT& comm):
	GlobalMatrixT(out, check_code),
	fComm(comm),
	fBuilder(NULL),
	fIsSymFactorized(false),
	fIsNumFactorized(false)
{
	const char caller[] = "SuperLU_DISTMatrixT::SuperLU_DISTMatrixT";

	fBuilder = new MSRBuilderT(false);
	if (!fBuilder) ExceptionT::OutOfMemory(caller);

	/* set up NULL structures */
	fA.Stype = SLU_NR_loc; /* distributed compressed row */
	fA.Dtype = SLU_D;  /* storing doubles */
	fA.Mtype = SLU_GE; /* general matrix */
	fA.nrow = 0;
	fA.ncol = 0;
	fA.Store = malloc(sizeof(NRformat_loc));
	if (!fA.Store) ExceptionT::OutOfMemory(caller);

	NRformat_loc *A = (NRformat_loc*) fA.Store;
	A->nnz_loc = 0;
	A->m_loc = 0;
	A->fst_row = 0;
	A->nzval = NULL;
	A->rowptr = NULL;
	A->colind = NULL;

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
#if __option (extended_errorcheck)
    foptions.PrintStat = YES;
#else
    foptions.PrintStat = NO;
#endif
	foptions.SymmetricMode = (symmetric) ? YES : NO;		
}

/* Destructor */	
SuperLU_DISTMatrixT::~SuperLU_DISTMatrixT(void)
{
	delete fBuilder;

	/* free the matrix */
	free(fA.Store);

	/* free upper and lower factors */
	if (fIsNumFactorized) {
		Destroy_SuperNode_Matrix(&fL);
		Destroy_CompCol_Matrix(&fU);
	}
}

/* add to structure */
void SuperLU_DISTMatrixT::AddEquationSet(const iArray2DT& eqnos) { fBuilder->AddGroup(eqnos); }
void SuperLU_DISTMatrixT::AddEquationSet(const RaggedArray2DT<int>& eqnos) { fBuilder->AddGroup(eqnos); }

/* set the internal matrix structure.
* NOTE: do not call Initialize() until equation topology has been set
* with AddEquationSet() for all equation sets */
void SuperLU_DISTMatrixT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	const char caller[] = "SuperLU_DISTMatrixT::Initialize";

	/* inherited */
	GlobalMatrixT::Initialize(tot_num_eq, loc_num_eq, start_eq);
	
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

#pragma message("what about these?")
	/* dimension work space */
	fperm_c.Dimension(fLocNumEQ);
	fperm_r.Dimension(fLocNumEQ);
	fetree.Dimension(fLocNumEQ);

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

	/* set update vector - global numbering */
	iArrayT activerows(fLocNumEQ);
	int n_update = fStartEQ; //OFFSET
	for (int i = 0; i < fLocNumEQ; i++)
		activerows[i] = n_update++;

	/* return the distributed SuperLU data structure */
	iArrayT rowptr;
	iArrayT colind;
	fBuilder->SetSuperLUData(activerows, iArrayT& rowptr, iArrayT& colind);
	frowptr = rowptr;
	fcolind = colind;
	fnzval.Dimension(colind.Length());

	/* set A matrix */
	fA.nrow = fTotNumEQ;
	fA.ncol = fTotNumEQ;
	NRformat_loc *A = (NRformat_loc*) fA.Store;
	A->nnz_loc = frowptr.Last();
	A->m_loc = fLocNumEQ;
	A->fst_row = fStartEQ - 1;
	A->nzval = fnzval.Pointer();
	A->rowptr = frowptr.Pointer();
	A->colind = fcolind.Pointer();

#pragma message("what about these?")
	/* scalings */
	fR.Dimension(fA.nrow);
	fC.Dimension(fA.ncol);

	/* output */
	fOut <<" Number of nonzeros in global matrix = "<< A->nnz <<"\n"<<endl;
	
	/* reset flags/options */
	fIsSymFactorized = false;
	fequed = 'N';
}

/* set all matrix values to 0.0 */
void SuperLU_DISTMatrixT::Clear(void)
{
	/* inherited */
	GlobalMatrixT::Clear();

	/* clear entries */
	dArrayT tmp;
	tmp.Alias(fnzval);
	tmp = 0.0;
	
	/* no equilibration */
	fequed = 'N';
	fIsNumFactorized = false;
}

/* add element group equations to the overall topology.
 * NOTE: assembly positions (equation numbers) = 1...fLocNumEQ
 * equations can be of fixed size (iArray2DT) or
 * variable length (RaggedArray2DT) */
void SuperLU_DISTMatrixT::AddEquationSet(const iArray2DT& eqset)
{
	fEqnos.AppendUnique (&eqset);
}

/* see AddEquationSet above */
void SuperLU_DISTMatrixT::AddEquationSet(const RaggedArray2DT<int>& eqset)
{
	fRaggedEqnos.AppendUnique (&eqset);
}

/* assemble the element contribution into the LHS matrix - assumes
* that elMat is square (n x n) and that eqnos is also length n.
* NOTE: assembly positions (equation numbers) = 1...fLocNumEQ */
void SuperLU_DISTMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos)
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

void SuperLU_DISTMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& row_eqnos,
	const ArrayT<int>& col_eqnos)
{
#pragma unused(elMat)
#pragma unused(row_eqnos)
#pragma unused(col_eqnos)
	ExceptionT::GeneralFail("SuperLU_DISTMatrixT::Assemble", "non-square not implemented");
}

void SuperLU_DISTMatrixT::Assemble(const nArrayT<double>& diagonal_elMat, const ArrayT<int>& eqnos)
{
#pragma unused(diagonal_elMat)
#pragma unused(eqnos)
	ExceptionT::GeneralFail("SuperLU_DISTMatrixT::Assemble", "diagonal not implemented");
}

/* number scope and reordering */
GlobalMatrixT::EquationNumberScopeT SuperLU_DISTMatrixT::EquationNumberScope(void) const { return kGlobal; }
bool SuperLU_DISTMatrixT::RenumberEquations(void) const { return false; }

GlobalMatrixT& SuperLU_DISTMatrixT::operator=(const GlobalMatrixT& rhs)
{
#pragma unused(rhs)
	ExceptionT::GeneralFail("SuperLU_DISTMatrixT::operator=", "not implemented");
	return *this;
}

/* return a clone of self */
GlobalMatrixT* SuperLU_DISTMatrixT::Clone(void) const
{
	/* not implemented */
	ExceptionT::GeneralFail("SuperLU_DISTMatrixT::operator=", "not implemented");
	return (GlobalMatrixT*) this;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* solution driver */
void SuperLU_DISTMatrixT::BackSubstitute(dArrayT& result)
{
	const char caller[] = "SuperLU_DISTMatrixT::BackSubstitute";
	
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
void SuperLU_DISTMatrixT::PrintAllPivots(void) const
{
// temp - not implemented yet. Maybe inappropriate
}

void SuperLU_DISTMatrixT::PrintZeroPivots(void) const
{
// temp - not implemented yet. Maybe inappropriate
}

void SuperLU_DISTMatrixT::PrintLHS(bool force) const
{
	if (!force || fCheckCode != GlobalMatrixT::kPrintLHS)
		return;

	fOut << "\nLHS matrix:\n\n";
	fOut << (*this) << "\n\n";
}

#endif /* __TAHOE_MPI__ */
#endif /* __SUPERLU_DIST__ */
