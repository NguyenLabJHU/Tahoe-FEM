/* $Id: AztecMatrixT.cpp,v 1.7.4.2 2002-10-20 18:07:46 paklein Exp $ */
/* created: paklein (08/10/1998) */

#include "AztecMatrixT.h"

/* library support options */
#ifdef __AZTEC__

#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>

#include "toolboxConstants.h"
#include "ExceptionT.h"

#include "Aztec_fe.h"
#include "iArray2DT.h"
#include "ElementMatrixT.h"
#include "RaggedArray2DT.h"

/* constructor */

using namespace Tahoe;

AztecMatrixT::AztecMatrixT(ifstreamT& in, ostream& out, int check_code):
	GlobalMatrixT(out, check_code),
	fInput(in)
{
	/* set and verify Aztec data structures */
	fAztec = new Aztec_fe(fInput, out);
	if (!fAztec) throw ExceptionT::kOutOfMemory;
}	

/* copy constructor */
AztecMatrixT::AztecMatrixT(const AztecMatrixT& source):
	GlobalMatrixT(source),
	fInput(source.fInput)
{
#pragma unused(source)
	cout << "\n AztecMatrixT::AztecMatrixT: not implemented" << endl;
	throw ExceptionT::kGeneralFail;
}

/* destuctor */
AztecMatrixT::~AztecMatrixT(void)
{
	/* free Aztec solver */
	delete fAztec;
	fAztec = NULL;
}

/* set matrix structure and allocate space.
* NOTE: do not call Initialize until all equations sets have been
* registered with SetStructure */
void AztecMatrixT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	/* inherited */
	GlobalMatrixT::Initialize(tot_num_eq, loc_num_eq, start_eq);

#ifndef __MPI__
	/* check */
	if (fTotNumEQ != fLocNumEQ)
	{
		cout << "\n AztecMatrixT::Initialize: no MPI: expecting the total number\n"
		     <<   "     of equations " << fTotNumEQ
		     << " to be equal to the local number of equations " << fLocNumEQ << endl;
		throw ExceptionT::kGeneralFail;
	}
#endif
	
	/* set-up Aztec matrix */
	fAztec->Initialize(loc_num_eq, start_eq);
	
	/* output statistics */
	int nonzerovals = fAztec->NumNonZeroValues();
	fOut << " Number of nonzero matrix values . . . . . . . . = ";
	fOut << nonzerovals << '\n';

	/* write Aztec options */
	fAztec->WriteAztecOptions(fOut);		
}

/* set all matrix values to 0.0 */
void AztecMatrixT::Clear(void)
{
	/* inherited */
	fAztec->Clear();
}

/* add element group equations to the overall topology.
* NOTE: assembly positions (equation numbers) = 1...fNumEQ
* equations can be of fixed size (iArray2DT) or
* variable length (RaggedArray2DT) */
void AztecMatrixT::AddEquationSet(const iArray2DT& eqset)
{
	/* inherited */
	fAztec->AddEquationSet(eqset);
	
	/* dimension workspace */
	int dim = eqset.MinorDim();
	if (dim > fValMat.Rows())
		fValMat.Dimension(dim);
}

void AztecMatrixT::AddEquationSet(const RaggedArray2DT<int>& eqset)
{
	/* inherited */
	fAztec->AddEquationSet(eqset);
	
	/* dimension workspace */
	int dim = eqset.MaxMinorDim();
	if (dim > fValMat.Rows())
		fValMat.Dimension(dim);
}

/* assemble the element contribution into the LHS matrix - assumes
* that elMat is square (n x n) and that eqnos is also length n.
*
* NOTE: assembly positions (equation numbers) = 1...fNumEQ */
void AztecMatrixT::Assemble(const ElementMatrixT& elMat, const nArrayT<int>& eqnos)
{
	/* element matrix format */
	ElementMatrixT::FormatT format = elMat.Format();

	if (format == ElementMatrixT::kDiagonal)
	{
		/* extract values for active equation numbers */
		fRowDexVec.Dimension(0);	
		fValVec.Dimension(0);
		int end_update = fStartEQ + fLocNumEQ - 1;
		for (int i = 0; i < eqnos.Length(); i++)
		{
			int eq = eqnos[i];
//			if (eq > 0)
			if (eq >= fStartEQ && eq <= end_update)
			{
				fRowDexVec.Append(eq - 1); //OFFSET
				fValVec.Append(elMat(i,i));
			}
		}
	
		/* assemble */
		int status;
		fAztec->AssembleDiagonals(fRowDexVec.Length(), fRowDexVec.Pointer(),
			fValVec.Pointer(), status);

		/* check completion */
		if (!status)
		{
			cout << "\n AztecMatrixT::Assemble: ERROR with equations:\n";
			cout << eqnos << endl;
			throw ExceptionT::kGeneralFail;
		}
	}
	else
	{
		/* equation numbers -> local active rows and all active columns */
		fRowDexVec.Dimension(0);
		fColDexVec.Dimension(0);
		int end_update = fStartEQ + fLocNumEQ - 1;
		for (int j = 0; j < eqnos.Length(); j++)
		{
			int eq = eqnos[j];
			if (eq >= fStartEQ && eq <= end_update)
				fRowDexVec.Append(j);
			if (eq > 0)
				fColDexVec.Append(j);
		}

		/* fill element matrix */
		if (elMat.Format() == ElementMatrixT::kSymmetricUpper)
			elMat.CopySymmetric();
	
		/* copy active block */
		int num_rows = fRowDexVec.Length();
		int num_cols = fColDexVec.Length();
		//TEMP - use nVariMatrixT
		dMatrixT activeblk(num_rows, num_cols, fValMat.Pointer());
		elMat.CopyBlock(fRowDexVec, fColDexVec, activeblk);
	
		/* active equation numbers -> global row numbers */
		for (int r = 0; r < num_rows; r++)
			fRowDexVec[r] = eqnos[fRowDexVec[r]] - 1; //OFFSET

		/* active equation numbers -> global col numbers */
		for (int c = 0; c < num_cols; c++)
			fColDexVec[c] = eqnos[fColDexVec[c]] - 1; //OFFSET

		/* row-by-row assembly */
		fValVec.Dimension(num_cols);
		int status = 1;
		for (int i = 0; i < num_rows && status; i++)
		{
			/* copy row values */
			activeblk.CopyRow(i, fValVec);
	
			/* assemble */
			fAztec->AssembleRow(fRowDexVec[i], num_cols, fColDexVec.Pointer(),
				fValVec.Pointer(), status);
		}
	
		/* check completion */
		if (!status)
		{
			cout << "\n AztecMatrixT::Assemble: ERROR with equations:\n";
			cout << eqnos << endl;
			throw ExceptionT::kGeneralFail;
		}
	}
}

void AztecMatrixT::Assemble(const ElementMatrixT& elMat, const nArrayT<int>& row_eqnos,
	const nArrayT<int>& col_eqnos)
{
	/* element matrix format */
	ElementMatrixT::FormatT format = elMat.Format();

	if (format == ElementMatrixT::kDiagonal)
	{
		cout << "\n AztecMatrixT::Assemble(m, r, c): cannot assemble diagonal matrix" << endl;
		throw ExceptionT::kGeneralFail;
	}
	else
	{
		int end_update = fStartEQ + fLocNumEQ - 1;

		/* equation numbers -> local active rows */
		fRowDexVec.Dimension(0);
		for (int j = 0; j < row_eqnos.Length(); j++)
		{
			int eq = row_eqnos[j];
			if (eq >= fStartEQ && eq <= end_update)
				fRowDexVec.Append(j);
		}

		/* equation numbers -> local active rows */
		fColDexVec.Dimension(0);
		for (int j = 0; j < col_eqnos.Length(); j++)
			if (col_eqnos[j] > 0)
				fColDexVec.Append(j);

		/* fill element matrix */
		if (elMat.Format() == ElementMatrixT::kSymmetricUpper)
			elMat.CopySymmetric();
	
		/* copy active block */
		int num_rows = fRowDexVec.Length();
		int num_cols = fColDexVec.Length();
		//TEMP - use nVariMatrixT
		dMatrixT activeblk(num_rows, num_cols, fValMat.Pointer());
		elMat.CopyBlock(fRowDexVec, fColDexVec, activeblk);
	
		/* active equation numbers -> global row numbers */
		for (int r = 0; r < num_rows; r++)
			fRowDexVec[r] = row_eqnos[fRowDexVec[r]] - 1; //OFFSET

		/* active equation numbers -> global col numbers */
		for (int c = 0; c < num_cols; c++)
			fColDexVec[c] = col_eqnos[fColDexVec[c]] - 1; //OFFSET

		/* row-by-row assembly */
		fValVec.Dimension(num_cols);
		int status = 1;
		for (int i = 0; i < num_rows && status; i++)
		{
			/* copy row values */
			activeblk.CopyRow(i, fValVec);
	
			/* assemble */
			fAztec->AssembleRow(fRowDexVec[i], num_cols, fColDexVec.Pointer(),
				fValVec.Pointer(), status);
		}
	
		/* check completion */
		if (!status)
		{
			cout << "\n AztecMatrixT::Assemble: ERROR with equations:\n";
			cout << " row:\n" << row_eqnos << '\n';
			cout << " col:\n" << col_eqnos << endl;
			throw ExceptionT::kGeneralFail;
		}
	}
}

/* number scope and reordering */
GlobalMatrixT::EquationNumberScopeT AztecMatrixT::EquationNumberScope(void) const
{
	return kGlobal;
}

bool AztecMatrixT::RenumberEquations(void) const { return false; }

/* assignment operator */
GlobalMatrixT& AztecMatrixT::operator=(const AztecMatrixT& rhs)
{
	cout <<  "\n AztecMatrixT::operator= : not implemented" << endl;
	throw ExceptionT::kGeneralFail;
	return *this;
}

/* assignment operator */
GlobalMatrixT& AztecMatrixT::operator=(const GlobalMatrixT& rhs)
{
#ifdef __NO_RTTI__
	cout << "\n AztecMatrixT::operator= : requires RTTI" << endl;
	throw ExceptionT::kGeneralFail;
#endif

	const AztecMatrixT* az = dynamic_cast<const AztecMatrixT*>(&rhs);
	if (!az) {
		cout << "\n AztecMatrixT::operator= : cast failed" << endl;
		throw ExceptionT::kGeneralFail;
	}
	return operator=(*az);
}
	
/* return a clone of self. Caller is responsible for disposing of the matrix */
GlobalMatrixT* AztecMatrixT::Clone(void) const
{
	AztecMatrixT* az = new AztecMatrixT(*this);
	return az;
}

/*************************************************************************
* Protected
*************************************************************************/

/* precondition matrix */
void AztecMatrixT::Factorize(void)
{
	/* preconditioning done during solve */
	fIsFactorized = 0; // undo GlobalMatrixT flag set
}
	
/* determine new search direction and put the results in result */
void AztecMatrixT::BackSubstitute(dArrayT& result)
{
	/* flag should not be set */
	if (fIsFactorized) throw ExceptionT::kGeneralFail;

	/* inherited - no initial guess */
	fAztec->Solve(result);
}

/* rank check functions */
void AztecMatrixT::PrintAllPivots(void) const
{
//not implemented
}

void AztecMatrixT::PrintZeroPivots(void) const
{
//not implemented
}

void AztecMatrixT::PrintLHS(void) const
{
	if (fCheckCode != GlobalMatrixT::kPrintLHS) return;

	/* inherited */
	fAztec->PrintNonZero(fOut);
}

/* library support options */
#endif /* __AZTEC__ */