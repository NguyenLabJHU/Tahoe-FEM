/* $Id: DiagonalMatrixT.cpp,v 1.20 2005-02-04 22:01:54 paklein Exp $ */
/* created: paklein (03/23/1997) */
#include "DiagonalMatrixT.h"
#include <iostream.h>
#include <iomanip.h>
#include "toolboxConstants.h"
#include "iArrayT.h"
#include "ElementMatrixT.h"

using namespace Tahoe;

/* constructor */
DiagonalMatrixT::DiagonalMatrixT(ostream& out, int check_code, AssemblyModeT mode):
	GlobalMatrixT(out, check_code),
	fIsFactorized(false)
{
	try { SetAssemblyMode(mode); }
	catch (ExceptionT::CodeT) { throw ExceptionT::kBadInputValue; }
}

/* copy constructor */
DiagonalMatrixT::DiagonalMatrixT(const DiagonalMatrixT& source):
	GlobalMatrixT(source),
	fMatrix(source.fMatrix),
	fMode(source.fMode),
	fIsFactorized(source.fIsFactorized)
{


}

/* set assemble mode */
void DiagonalMatrixT::SetAssemblyMode(AssemblyModeT mode)
{
	/* set */
	fMode = mode;
	
	/* check */
	if (fMode != kNoAssembly &&
	    fMode != kDiagOnly   &&
	    fMode != kAbsRowSum) throw ExceptionT::kGeneralFail;
}

/* set the internal matrix structure.
* NOTE: do not call Initialize() equation topology has been set
* with AddEquationSet() for all equation sets */
void DiagonalMatrixT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	/* inherited */
	GlobalMatrixT::Initialize(tot_num_eq, loc_num_eq, start_eq);

	/* allocate work space */
	fMatrix.Dimension(fLocNumEQ);
	fIsFactorized = false;
}

/* set all matrix values to 0.0 */
void DiagonalMatrixT::Clear(void)
{
	/* inherited */
	GlobalMatrixT::Clear();

	/* clear data */
	fMatrix = 0.0;
	fIsFactorized = false;
}

/* add element group equations to the overall topology.
* NOTE: assembly positions (equation numbers) = 1...fNumEQ
* equations can be of fixed size (iArray2DT) or
* variable length (RaggedArray2DT) */
void DiagonalMatrixT::AddEquationSet(const iArray2DT& eqset)
{
#pragma unused(eqset)
// no equation data needed
}

void DiagonalMatrixT::AddEquationSet(const RaggedArray2DT<int>& eqset)
{
#pragma unused(eqset)
// no equation data needed
}

/* assemble the element contribution into the LHS matrix - assumes
* that elMat is square (n x n) and that eqnos is also length n.
* NOTE: assembly positions (equation numbers) = 1...fNumEQ */
void DiagonalMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos)
{
	if (elMat.Format() == ElementMatrixT::kDiagonal)
	{
		/* from diagonal only */
		const double* pelMat = elMat.Pointer();
		int inc = elMat.Rows() + 1;

		int numvals = eqnos.Length();
		for (int i = 0; i < numvals; i++)
		{
			int eq = eqnos[i];
		
			/* active dof */
			if (eq-- > 0) fMatrix[eq] += *pelMat;
	
			pelMat += inc;
		}	
	}
	else
	{
		/* assemble modes */
		switch (fMode)
		{
			case kDiagOnly:
			{		
				/* assemble only the diagonal values */
				for (int i = 0; i < eqnos.Length(); i++)
					if (eqnos[i] > 0)
						fMatrix[eqnos[i] - 1] += elMat(i,i);	
				break;
			}
			case kAbsRowSum:
			{
				/* copy to full symmetric */
				if (elMat.Format() == ElementMatrixT::kSymmetricUpper)
					elMat.CopySymmetric();

				/* assemble sum of row absolute values */
				int numrows = eqnos.Length();
				for (int i = 0; i < numrows; i++)
					if (eqnos[i] > 0)
					{
						const double* prow = elMat.Pointer(i);
						double sum = 0.0;					
						for (int j = 0; j < numrows; j++)
						{				
							sum  += fabs(*prow);
							prow += numrows;
						}

						fMatrix[eqnos[i] - 1] += sum;
					}
			
				break;
			}
			default:			
				cout << "\n DiagonalMatrixT::Assemble: cannot assemble mode: " << fMode << endl;
				throw ExceptionT::kGeneralFail;
		}
	}
}

void DiagonalMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& row_eqnos,
	const ArrayT<int>& col_eqnos)
{
	/* pick out diagonal values */
	for (int row = 0; row < row_eqnos.Length(); row++)
		for (int col = 0; col < col_eqnos.Length(); col++)
			if (row_eqnos[row] == col_eqnos[col])
			{
				int eqno = row_eqnos[row] - 1;
				if (eqno > -1)
					fMatrix[eqno] += elMat(row, col);
			}
}

void DiagonalMatrixT::Assemble(const nArrayT<double>& diagonal_elMat, const ArrayT<int>& eqnos)
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (diagonal_elMat.Length() != eqnos.Length()) throw ExceptionT::kSizeMismatch;
#endif

	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eqno = eqnos[i] - 1;
		if (eqno > -1)
			fMatrix[eqno] += diagonal_elMat[i];
	}
}

/* fetch values */
void DiagonalMatrixT::DisassembleDiagonal(dArrayT& diagonals, 
	const nArrayT<int>& eqnos) const
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (diagonals.Length() != eqnos.Length()) throw ExceptionT::kSizeMismatch;
#endif

	for (int i = 0; i < eqnos.Length(); i++)
	{
		/* ignore requests for inactive equations */	
		int eq = eqnos[i];	
		if (eq-- > 0)
			diagonals[i] = fMatrix[eq];
		else
			diagonals[i] = 0.0;
	}
}

/* number scope and reordering */
GlobalMatrixT::EquationNumberScopeT DiagonalMatrixT::EquationNumberScope(void) const
{
	return kLocal;
}

bool DiagonalMatrixT::RenumberEquations(void) const { return false; }

/* assignment operator */
DiagonalMatrixT& DiagonalMatrixT::operator=(const DiagonalMatrixT& rhs)
{
	/* no copies of self */
	if (this == &rhs) return *this;

	/* inherited */
	GlobalMatrixT::operator=(rhs);

	fMatrix = rhs.fMatrix;
	fMode   = rhs.fMode;
	fIsFactorized = rhs.fIsFactorized;

	return *this;
}

/** return a clone of self */
GlobalMatrixT* DiagonalMatrixT::Clone(void) const
{
	DiagonalMatrixT* new_mat = new DiagonalMatrixT(*this);
	return new_mat;
}

/* matrix-vector product */
void DiagonalMatrixT::Multx(const dArrayT& x, dArrayT& b) const
{
	b = x;
	if (fIsFactorized)
		b /= fMatrix;
	else
		b *= fMatrix;
}

/* vector-matrix-vector product */
double DiagonalMatrixT::MultmBn(const dArrayT& m, const dArrayT& n) const
{
#if __option(extended_errorcheck)
	if (m.Length() != fMatrix.Length() || n.Length() != fMatrix.Length())
		ExceptionT::SizeMismatch("DiagonalMatrixT::MultmBn");
#endif

	const double* pm = m.Pointer();
	const double* pn = n.Pointer();
	const double* pM = fMatrix.Pointer();
	int length = fMatrix.Length();
	double mBn = 0.0;
	for (int i = 0; i < length; i++)
		mBn += (*pm++)*(*pM++)*(*pn++);
	return mBn;
}

/**************************************************************************
* Protected
**************************************************************************/

/* precondition matrix */
void DiagonalMatrixT::Factorize(void)
{
	/* quick exit */
	if (fIsFactorized)
		return;
	else {
	
		double* pMatrix = fMatrix.Pointer();

		/* inverse of a diagonal matrix */
		bool first = true;
		for (int i = 0; i < fLocNumEQ; i++)
		{
			/* check for zero pivots */
			if (fabs(*pMatrix) > kSmall)
				*pMatrix = 1.0/(*pMatrix);
			else
			{
				if (first)
				{
					cout << "\n DiagonalMatrixT::Factorize: WARNING: skipping small pivots. See out file."
					     << endl;
				
					fOut << "\n DiagonalMatrixT::Factorize: small pivots\n";
					fOut << setw(kIntWidth) << "eqn"
					     << setw(OutputWidth(fOut, pMatrix)) << "pivot" << '\n';
					first = false;		
				}
			
				fOut << setw(kIntWidth) << i+1
				     << setw(OutputWidth(fOut, pMatrix)) << *pMatrix << '\n';
			}
	
			/* next */
			pMatrix++;
		}

		/* flush */
		if (!first) fOut.flush();
		
		/* set flag */
		fIsFactorized = true;
	}
}
	
/* solution driver */
void DiagonalMatrixT::BackSubstitute(dArrayT& result)
{
	/* checks */
	if (result.Length() != fLocNumEQ || !fIsFactorized) throw ExceptionT::kGeneralFail;
	
	double* presult = result.Pointer();
	double* pMatrix = fMatrix.Pointer();
	
	for (int i = 0; i < fLocNumEQ; i++)
		*presult++ *= *pMatrix++;
}

/* check functions */
void DiagonalMatrixT::PrintAllPivots(void) const
{
	if (fCheckCode != GlobalMatrixT::kAllPivots)
		return;

	fOut << "\nAll pivots:\n\n";
	fOut << fMatrix << "\n\n";
}

void DiagonalMatrixT::PrintZeroPivots(void) const
{
	if (fCheckCode != GlobalMatrixT::kZeroPivots) return;
	int d_width = OutputWidth(fOut, fMatrix.Pointer());

	int firstline = 1;
	for (int i = 0; i < fLocNumEQ; i++)
	{
		double pivot = fMatrix[i];
		
		if (pivot < kSmall)
		{
			if (firstline)
			{
				fOut << "\nZero or negative pivots:\n\n";
				firstline = 0;
			}
		
			fOut << setw(kIntWidth) << i + 1;
			fOut << setw(d_width) << pivot << '\n';
		}
	}
	
	if (!firstline) fOut << '\n';
}

void DiagonalMatrixT::PrintLHS(bool force) const
{
	if (!force && fCheckCode != GlobalMatrixT::kPrintLHS) return;
		
	fOut << "\nLHS matrix:\n\n";
	fOut << fMatrix << "\n\n";
}
