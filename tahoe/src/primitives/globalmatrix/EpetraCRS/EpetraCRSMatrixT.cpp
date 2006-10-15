/* $Id: EpetraCRSMatrixT.cpp,v 1.1 2006-10-15 04:14:14 paklein Exp $ */
#include "EpetraCRSMatrixT.h"

/* library support options */
#ifdef __TRILINOS__

#include "MSRBuilderT.h"
#include "dMatrixT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "ElementMatrixT.h"
#include "CommunicatorT.h"

/* Epetra headers */
#include "Epetra_ConfigDefs.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#ifdef __TAHOE_MPI__
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

using namespace Tahoe;

/* constructor */
EpetraCRSMatrixT::EpetraCRSMatrixT(ostream& out, int check_code, const CommunicatorT& comm):
	GlobalMatrixT(out, check_code, comm),
	fBuilder(NULL),
	fIsSymFactorized(false),
	fIsNumFactorized(false)
{
	const char caller[] = "EpetraCRSMatrixT::EpetraCRSMatrixT";

	fBuilder = new MSRBuilderT(false);
	if (!fBuilder) ExceptionT::OutOfMemory(caller);
}

/* copy constructor */
EpetraCRSMatrixT::EpetraCRSMatrixT(const EpetraCRSMatrixT& rhs):
	GlobalMatrixT(rhs),
	fBuilder(NULL)
{
	EpetraCRSMatrixT::operator=(rhs);
}

/* Destructor */	
EpetraCRSMatrixT::~EpetraCRSMatrixT(void)
{
	delete fBuilder;
}

/* add to structure */
void EpetraCRSMatrixT::AddEquationSet(const iArray2DT& eqnos) { fBuilder->AddGroup(eqnos); }
void EpetraCRSMatrixT::AddEquationSet(const RaggedArray2DT<int>& eqnos) { fBuilder->AddGroup(eqnos); }

/* set the internal matrix structure.
* NOTE: do not call Initialize() until equation topology has been set
* with AddEquationSet() for all equation sets */
void EpetraCRSMatrixT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	const char caller[] = "EpetraCRSMatrixT::Initialize";

	/* inherited */
	GlobalMatrixT::Initialize(tot_num_eq, loc_num_eq, start_eq);

	/* set update vector - global numbering */
	factive.Dimension(fLocNumEQ);
	int n_update = fStartEQ; //OFFSET
	for (int i = 0; i < fLocNumEQ; i++)
		factive[i] = n_update++;

	/* return the distributed structure data */
	iArrayT active; active.Alias(factive);
	iArrayT rowptr;
	iArrayT colind;
	fBuilder->SetSuperLUData(active, rowptr, colind);
	frowptr = rowptr;
	fcolind = colind;
	fnzval.Dimension(colind.Length());

	/* reset flags/options */
	fIsSymFactorized = false;
	fIsNumFactorized = false;
}

/* set all matrix values to 0.0 */
void EpetraCRSMatrixT::Clear(void)
{
	/* inherited */
	GlobalMatrixT::Clear();

	/* clear entries */
	dArrayT tmp;
	tmp.Alias(fnzval);
	tmp = 0.0;
	
	/* no equilibration */
	fIsNumFactorized = false;
}

/* assemble the element contribution into the LHS matrix - assumes
* that elMat is square (n x n) and that eqnos is also length n.
* NOTE: assembly positions (equation numbers) = 1...fLocNumEQ */
void EpetraCRSMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos)
{
	const char caller[] = "EpetraCRSMatrixT::Assemble";

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

void EpetraCRSMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& row_eqnos,
	const ArrayT<int>& col_eqnos)
{
#pragma unused(elMat)
#pragma unused(row_eqnos)
#pragma unused(col_eqnos)
	ExceptionT::GeneralFail("EpetraCRSMatrixT::Assemble", "non-square not implemented");
}

void EpetraCRSMatrixT::Assemble(const nArrayT<double>& diagonal_elMat, const ArrayT<int>& eqnos)
{
#pragma unused(diagonal_elMat)
#pragma unused(eqnos)
	ExceptionT::GeneralFail("EpetraCRSMatrixT::Assemble", "diagonal not implemented");
}

/* number scope and reordering */
GlobalMatrixT::EquationNumberScopeT EpetraCRSMatrixT::EquationNumberScope(void) const { return kGlobal; }
bool EpetraCRSMatrixT::RenumberEquations(void) const { return false; }

EpetraCRSMatrixT& EpetraCRSMatrixT::operator=(const EpetraCRSMatrixT&)
{
	ExceptionT::GeneralFail("EpetraCRSMatrixT::operator=", "not implemented");
	return *this;
}

/* return a clone of self */
GlobalMatrixT* EpetraCRSMatrixT::Clone(void) const {
	return new EpetraCRSMatrixT(*this);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* solution driver */
void EpetraCRSMatrixT::BackSubstitute(dArrayT& result)
{
	const char caller[] = "EpetraCRSMatrixT::BackSubstitute";

	/* set-up Epetra structures */
#ifdef __TAHOE_MPI__
	Epetra_MpiComm epetra_comm(fComm);
#else
	Epetra_SerialComm epetra_comm;
#endif
	Epetra_Map epetra_map(fTotNumEQ, fLocNumEQ, factive.Pointer(), 1, epetra_comm);
	iArrayT row_count;
	// fill values from frowptr
	//
	//
	Epetra_CrsMatrix epetra_matrix(Copy, epetra_map, row_count.Pointer(), true);

	/* copy data into epetra_matrix */
	for (int i = 0; i < factive.Length(); i++) {
		int offset = frowptr[i];
		epetra_matrix.InsertGlobalValues(factive[i], row_count[i], fnzval.Pointer(offset), fcolind.Pointer(offset));
		
		// fcolind - 0 or 1 numbering?

	}  
	epetra_matrix.FillComplete();

	/* call solver */

	/* always fully factorized on exit */
	fIsSymFactorized = true;
	fIsNumFactorized = true;
}

/* check functions */
void EpetraCRSMatrixT::PrintAllPivots(void) const
{
//not implemented
}

void EpetraCRSMatrixT::PrintZeroPivots(void) const
{
//not implemented
}

void EpetraCRSMatrixT::PrintLHS(bool force) const
{
	if (!force && fCheckCode != GlobalMatrixT::kPrintLHS) return;

#if 0
	const char caller[] = "EpetraCRSMatrixT::PrintLHS";
	fOut << "\n " << caller << '\n';

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
#endif
}

#endif /* __TRILINOS__ */
