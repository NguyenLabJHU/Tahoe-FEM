/* $Id: SuperLU_DISTMatrixT.cpp,v 1.8 2005-08-01 03:26:30 paklein Exp $ */
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
#include "CommunicatorT.h"

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

using namespace Tahoe;

/* constructor */
SuperLU_DISTMatrixT::SuperLU_DISTMatrixT(ostream& out, int check_code, bool print_stat, 
	IterRefine_t refine, const CommunicatorT& comm):
	GlobalMatrixT(out, check_code, comm),
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

	/* create blank data structures */
	fScalePermstruct.DiagScale = NOEQUIL;
	fScalePermstruct.R = NULL;
	fScalePermstruct.C = NULL;
	fScalePermstruct.perm_r = NULL;
	fScalePermstruct.perm_c = NULL;
	
	fLUstruct.etree = NULL;
	fLUstruct.Glu_persist = NULL;
	fLUstruct.Llu = NULL;	

    /* Set the default input options:
        options.Fact = DOFACT;
        options.Equil = YES;
        options.ColPerm = MMD_AT_PLUS_A;
        options.RowPerm = LargeDiag;
        options.ReplaceTinyPivot = YES;
        options.Trans = NOTRANS;
        options.IterRefine = DOUBLE;
        options.SolveInitialized = NO;
        options.RefineInitialized = NO;
        options.PrintStat = YES;
     */
    set_default_options_dist(&foptions);
    foptions.PrintStat = (print_stat) ? YES : NO;
    foptions.IterRefine = refine;

	/* initialize the process grid - divide as square */
	int size = fComm.Size();
	int nprow = int(sqrt(double(size)));
	nprow = (nprow < 1) ? 1 : nprow;
	int npcol = size/nprow;
	if (npcol*nprow != size)
		ExceptionT::GeneralFail(caller, "number of processes must be even %d", size);
    superlu_gridinit(fComm.Comm(), nprow, npcol, &fgrid);
}

/* copy constructor */
SuperLU_DISTMatrixT::SuperLU_DISTMatrixT(const SuperLU_DISTMatrixT& rhs):
	GlobalMatrixT(rhs),
	fBuilder(NULL)
{
	SuperLU_DISTMatrixT::operator=(rhs);
}

/* Destructor */	
SuperLU_DISTMatrixT::~SuperLU_DISTMatrixT(void)
{
	delete fBuilder;

	/* free the matrix struct */
	SUPERLU_FREE(fA.Store);

	/* free upper and lower factors */
	FreeScalePermstruct(fScalePermstruct);
	FreeLUstruct(fLocNumEQ, fgrid, fLUstruct);

	/* free information about solution phase */
    if (foptions.SolveInitialized)
        dSolveFinalize(&foptions, &fSOLVEstruct);

	/* release process grid */
    superlu_gridexit(&fgrid);
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

	/* free existing memory */
	FreeLUstruct(fLocNumEQ, fgrid, fLUstruct);
	FreeScalePermstruct(fScalePermstruct);

	/* free information about solution phase */
    if (foptions.SolveInitialized) {
        dSolveFinalize(&foptions, &fSOLVEstruct);
        foptions.SolveInitialized = NO;
    }

	/* inherited */
	GlobalMatrixT::Initialize(tot_num_eq, loc_num_eq, start_eq);

	/* redimension */
    ScalePermstructInit(fTotNumEQ, fTotNumEQ, &fScalePermstruct);
    LUstructInit(fTotNumEQ, fTotNumEQ, &fLUstruct);

	/* set update vector - global numbering */
	iArrayT activerows(fLocNumEQ);
	int n_update = fStartEQ; //OFFSET
	for (int i = 0; i < fLocNumEQ; i++)
		activerows[i] = n_update++;

	/* return the distributed SuperLU data structure */
	iArrayT rowptr;
	iArrayT colind;
	fBuilder->SetSuperLUData(activerows, rowptr, colind);
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

	/* reset flags/options */
	fIsSymFactorized = false;
	fIsNumFactorized = false;
}

/* write information to output stream after SuperLU_DISTMatrixT::Initialize
 * has been called */
void SuperLU_DISTMatrixT::Info(ostream& out)
{
	/* inherited */
	GlobalMatrixT::Info(out);

	/* output */
	NRformat_loc *A = (NRformat_loc*) fA.Store;
	out <<" Number of nonzeros in local global matrix = "<< A->nnz_loc <<"\n"<<endl;	
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
	fIsNumFactorized = false;

	/* free information about solution phase */
    if (foptions.SolveInitialized) {
        dSolveFinalize(&foptions, &fSOLVEstruct);
        foptions.SolveInitialized = NO;
    }
}

/* assemble the element contribution into the LHS matrix - assumes
* that elMat is square (n x n) and that eqnos is also length n.
* NOTE: assembly positions (equation numbers) = 1...fLocNumEQ */
void SuperLU_DISTMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos)
{
	const char caller[] = "SuperLU_DISTMatrixT::Assemble";

	/* element matrix format */
	ElementMatrixT::FormatT format = elMat.Format();

	/* two cases: element matrix is diagonal, or it's not. */
	int end_update = fStartEQ + fLocNumEQ - 1;
	if (format == ElementMatrixT::kDiagonal)
	{
		/* diagonal entries only */
		const double *pelMat = elMat.Pointer();
		int inc = elMat.Rows() + 1; /* offset between diag entries are */
		int nee = eqnos.Length();
		for (int i = 0; i < nee; ++i) {
			int eq = eqnos[i];
			if (eq >= fStartEQ && eq <= end_update) /* active eqn */ {
				eq--;
				double* a = (*this)(eq,eq);
				if (a)
					*a += *pelMat;
				else
					ExceptionT::OutOfRange(caller);
			}
			pelMat += inc;
		}
	}
	else if (format == ElementMatrixT::kNonSymmetric || 
             format == ElementMatrixT::kSymmetric ||
             format == ElementMatrixT::kSymmetricUpper )
	{
		/* fill matrix */
		if (format != ElementMatrixT::kNonSymmetric)
			elMat.CopySymmetric();

		int nee = eqnos.Length();  // number of equations for element
		for (int col = 0; col < nee; ++col)
		{
			int ceqno = eqnos[col] - 1;
			if (ceqno > -1) /* active eqn */ {
				for (int row = 0; row < nee; ++row) {
					int reqno = eqnos[row];
					if (reqno >= fStartEQ && reqno <= end_update) /* active eqn */ {
						reqno--;
						double* a = (*this)(reqno,ceqno);
						if (a)
							*a += 0.5*(elMat(row,col) + elMat(col,row));
						else {
							iArrayT tmp;
							tmp.Alias(eqnos);
							fOut << "\n " << caller << ": bad eqnos = " << tmp.no_wrap() << endl;
							ExceptionT::OutOfRange(caller);
						}
					}
				}
			}
		}
	}
	else
		ExceptionT::GeneralFail(caller, "unsupported element matrix format %d", format);
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

SuperLU_DISTMatrixT& SuperLU_DISTMatrixT::operator=(const SuperLU_DISTMatrixT&)
{
	ExceptionT::GeneralFail("SuperLU_DISTMatrixT::operator=", "not implemented");
	return *this;
}

/* return a clone of self */
GlobalMatrixT* SuperLU_DISTMatrixT::Clone(void) const {
	return new SuperLU_DISTMatrixT(*this);
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

    /* Initialize the statistics variables. */
	SuperLUStat_t stat;    
    PStatInit(&stat);

	/* make copy since solver may permute */
	fcolind2 = fcolind;
	frowptr2 = frowptr;

	/* call SuperLU */
    int info;
    int nrhs = 1;
	double berr;
    pdgssvx(&foptions, &fA, &fScalePermstruct, result.Pointer(), result.Length(), nrhs, &fgrid,
		     &fLUstruct, &fSOLVEstruct, &berr, &stat, &info);

	/* restore */
	fcolind = fcolind2;
	frowptr = frowptr2;

	/* check results */
	if (info != 0)
		ExceptionT::BadJacobianDet(caller, "pdgssvx_ABglobal return %d with backward componentwise error %g", info, berr);

	/* report statistics */
    if (foptions.PrintStat) PStatPrint(&foptions, &stat, &fgrid);
    PStatFree(&stat);

	/* always fully factorized on exit */
	foptions.SolveInitialized = yes_no_t(YES);
	fIsSymFactorized = true;
	fIsNumFactorized = true;
}

/* check functions */
void SuperLU_DISTMatrixT::PrintAllPivots(void) const
{
//not implemented
}

void SuperLU_DISTMatrixT::PrintZeroPivots(void) const
{
//not implemented
}

void SuperLU_DISTMatrixT::PrintLHS(bool force) const
{
	if (!force && fCheckCode != GlobalMatrixT::kPrintLHS) return;

	const char caller[] = "SuperLU_DISTMatrixT::PrintLHS";
	fOut << "\n " << caller << '\n';
	
	NRformat_loc *A = (NRformat_loc*) fA.Store;
	fOut << "number of nonzero = " << A->nnz_loc << '\n';
	fOut << "rows on this processor = " << A->m_loc << '\n';
	fOut << "first row on this processor = " << A->fst_row << '\n';
	
	/* not allocated */
	if (A->nzval == NULL || A->rowptr == NULL || A->colind == NULL) {
		fOut << endl;
		ExceptionT::GeneralFail(caller, "storage is NULL");
	}

	/* sparse matrix format */
	iArrayT tmp;
	tmp.Alias(frowptr);
	fOut << "row pointers:\n" << tmp.wrap(10) << '\n';
	fOut << "col indicies:\n";
	for (int i = 0; i < A->m_loc; i++) {
		int length = frowptr[i+1] - frowptr[i];
		tmp.Alias(length, fcolind.Pointer(frowptr[i]));
		fOut << tmp.no_wrap() << '\n';
	}

	fOut << "LHS: {r, c, v}: \n";
	int dim = A->m_loc;
	const double* nzval_ = (const double*) A->nzval;
	const int* rowptr_ = A->rowptr;
	const int* colind_ = A->colind;
	for (int i = 0; i < dim; i++) {
	
		int index = rowptr_[i];
		int count = rowptr_[i+1] - index;

		const int* col = colind_ + index;
		const double* val = nzval_ + index;
	
		for (int j = 0; j < count; j++)
			fOut << i+fStartEQ << " " << (*col++)+1 << " " << *val++ << '\n';
	}
	fOut << endl;
}

void SuperLU_DISTMatrixT::FreeLUstruct(int dim, gridinfo_t& grid, LUstruct_t& lu_struct) const
{
	/* only free if allocated */
	if (lu_struct.Glu_persist) {
		Destroy_LU(dim, &grid, &lu_struct);
		LUstructFree(&lu_struct);
	}
	
	/* clear all pointers */
	lu_struct.etree = NULL;
	lu_struct.Glu_persist = NULL;
	lu_struct.Llu = NULL;
}

void SuperLU_DISTMatrixT::FreeScalePermstruct(ScalePermstruct_t& scale_struct) const
{
	/* if allocated */
	if (scale_struct.perm_r)
		ScalePermstructFree(&scale_struct);

	/* clear all pointers */
	scale_struct.R = NULL;
	scale_struct.C = NULL;
	scale_struct.perm_r = NULL;
	scale_struct.perm_c = NULL;
}

#endif /* __TAHOE_MPI__ */
#endif /* __SUPERLU_DIST__ */
