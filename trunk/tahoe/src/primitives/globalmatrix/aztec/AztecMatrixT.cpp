/* $Id: AztecMatrixT.cpp,v 1.24 2005-04-13 21:50:27 paklein Exp $ */
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

using namespace Tahoe;

/* constructor */
AztecMatrixT::AztecMatrixT(ostream& out, int check_code, const CommunicatorT& comm, 
	const ParameterListT& parameters):
	GlobalMatrixT(out, check_code, comm)
{
	/* set and verify Aztec data structures */
	fAztec = new Aztec_fe(parameters, out, comm);
	if (!fAztec) throw ExceptionT::kOutOfMemory;
}	

/* copy constructor */
AztecMatrixT::AztecMatrixT(const AztecMatrixT& source):
	GlobalMatrixT(source)
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

#ifndef __TAHOE_MPI__
	/* check */
	if (fTotNumEQ != fLocNumEQ)
		ExceptionT::GeneralFail("AztecMatrixT::Initialize",
			"no MPI: total equations %d != local equations %d", fTotNumEQ, fLocNumEQ);
#endif
	
	/* set-up Aztec matrix */
	fAztec->Initialize(loc_num_eq, start_eq);
}

/* write information to output stream after GlobalMatrixT::Initialize
 * has been called */
void AztecMatrixT::Info(ostream& out)
{
	/* inherited */
	GlobalMatrixT::Info(out);
	
	/* output statistics */
	int nonzerovals = fAztec->NumNonZeroValues();
	out << " Number of nonzero matrix values . . . . . . . . = ";
	out << nonzerovals << '\n';

	/* write Aztec options */
	fAztec->WriteAztecOptions(out);
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
void AztecMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos)
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
			iArrayT tmp;
			tmp.Alias(eqnos);
			cout << "\n AztecMatrixT::Assemble: ERROR with equations:\n";
			cout << tmp << endl;
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
			iArrayT tmp;
			tmp.Alias(eqnos);
			cout << "\n AztecMatrixT::Assemble: ERROR with equations:\n";
			cout << tmp << endl;
			throw ExceptionT::kGeneralFail;
		}
	}
}

void AztecMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& row_eqnos,
	const ArrayT<int>& col_eqnos)
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
			iArrayT tmp;
			cout << "\n AztecMatrixT::Assemble: ERROR with equations:\n";
			tmp.Alias(row_eqnos);
			cout << " row:\n" << tmp << '\n';
			tmp.Alias(col_eqnos);
			cout << " col:\n" << tmp << endl;
			throw ExceptionT::kGeneralFail;
		}
	}
}

void AztecMatrixT::Assemble(const nArrayT<double>& diagonal_elMat, const ArrayT<int>& eqnos)
{
#pragma unused(diagonal_elMat)
#pragma unused(eqnos)

	ExceptionT::GeneralFail("AztecMatrixT::Assemble", "not implemented");
}

/* number scope and reordering */
GlobalMatrixT::EquationNumberScopeT AztecMatrixT::EquationNumberScope(void) const
{
	return kGlobal;
}

bool AztecMatrixT::RenumberEquations(void) const { return false; }

/* assignment operator */
AztecMatrixT& AztecMatrixT::operator=(const AztecMatrixT&)
{
	const char caller[] = "AztecMatrixT::operator=";
	ExceptionT::GeneralFail(caller, "not implemented");	
	return *this;
}
	
/* return a clone of self. Caller is responsible for disposing of the matrix */
GlobalMatrixT* AztecMatrixT::Clone(void) const {
	return new AztecMatrixT(*this);
}

/*************************************************************************
 * Protected
 *************************************************************************/
	
/* determine new search direction and put the results in result */
void AztecMatrixT::BackSubstitute(dArrayT& result)
{
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

void AztecMatrixT::PrintLHS(bool force) const
{
	if (!force && fCheckCode != GlobalMatrixT::kPrintLHS) return;

	/* inherited */
	fAztec->PrintNonZero(fOut);
}

/* library support options */
#endif /* __AZTEC__ */
