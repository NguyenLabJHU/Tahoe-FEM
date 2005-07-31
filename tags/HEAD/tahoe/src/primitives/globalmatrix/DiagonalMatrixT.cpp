/* $Id: DiagonalMatrixT.cpp,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (03/23/1997)                                          */
/* Virtual base class for all global matrix objects                       */

#include "DiagonalMatrixT.h"
#include <iostream.h>
#include <iomanip.h>
#include "Constants.h"
#include "iArrayT.h"
#include "ElementMatrixT.h"

/* constructor */
DiagonalMatrixT::DiagonalMatrixT(ostream& out, int check_code, AssemblyModeT mode):
	GlobalMatrixT(out, check_code)
{
	try { SetAssemblyMode(mode); }
	catch (int) { throw eBadInputValue; }
}

/* set assemble mode */
void DiagonalMatrixT::SetAssemblyMode(AssemblyModeT mode)
{
	/* set */
	fMode = mode;
	
	/* check */
	if (fMode != kNoAssembly &&
	    fMode != kDiagOnly   &&
	    fMode != kAbsRowSum) throw eGeneralFail;
}

/* set the internal matrix structure.
* NOTE: do not call Initialize() equation topology has been set
* with AddEquationSet() for all equation sets */
void DiagonalMatrixT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	/* inherited */
	GlobalMatrixT::Initialize(tot_num_eq, loc_num_eq, start_eq);

	/* allocate work space */
	fMatrix.Allocate(fLocNumEQ);
}

/* set all matrix values to 0.0 */
void DiagonalMatrixT::Clear(void)
{
	/* inherited */
	GlobalMatrixT::Clear();

	/* clear data */
	fMatrix = 0.0;
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
void DiagonalMatrixT::Assemble(const ElementMatrixT& elMat, const iArrayT& eqnos)
{
	if (elMat.Format() == ElementMatrixT::kDiagonal)
	{
		/* from diagonal only */
		double* pelMat = elMat.Pointer();
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
						double* prow = elMat.Pointer(i);
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
				/* no assembly mode specified */
				throw eGeneralFail;
		}
	}
}

/* assignment operator */
GlobalMatrixT& DiagonalMatrixT::operator=(const GlobalMatrixT& RHS)
{
#ifdef __NO_RTTI__
	const DiagonalMatrixT* pRHS = (const DiagonalMatrixT*) (&RHS);
#else
	const DiagonalMatrixT* pRHS = dynamic_cast<const DiagonalMatrixT*>(&RHS);
	if (!pRHS) throw eGeneralFail;
#endif

	fMatrix = pRHS->fMatrix;

	/* inherited */
	return GlobalMatrixT::operator=(RHS);
}

/* number scope and reordering */
GlobalMatrixT::EquationNumberScopeT DiagonalMatrixT::EquationNumberScope(void) const
{
	return kLocal;
}

bool DiagonalMatrixT::RenumberEquations(void) const { return false; }

/**************************************************************************
* Protected
**************************************************************************/

/* precondition matrix */
void DiagonalMatrixT::Factorize(void)
{
	double* pMatrix = fMatrix.Pointer();

	/* inverse of a diagonal matrix */
	bool first = true;
	for (int i = 0; i < fLocNumEQ; i++)
	{
		/* check zero pivots */
		if (fabs(*pMatrix) < kSmall)
		{
			if (first)
			{
				cout << "\n DiagonalMatrixT::Factorize: WARNING: small pivots. See out file."
				     << endl;
				
				fOut << "\n DiagonalMatrixT::Factorize: small pivots\n";
				fOut << setw(kIntWidth) << "eqn"
				     << setw(OutputWidth(fOut, pMatrix)) << "pivot" << '\n';
				first = false;		
			}
		
			fOut << setw(kIntWidth) << i+1
			     << setw(OutputWidth(fOut, pMatrix)) << *pMatrix << '\n';
		}
	
		*pMatrix = 1.0/(*pMatrix);
		pMatrix++;
	}

	/* flush */
	if (!first) fOut.flush();

	fIsFactorized = 1;
}
	
/* solution driver */
void DiagonalMatrixT::BackSubstitute(dArrayT& result)
{
	/* checks */
	if (result.Length() != fLocNumEQ || !fIsFactorized) throw eGeneralFail;
	
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

void DiagonalMatrixT::PrintLHS(void) const
{
	if (fCheckCode != GlobalMatrixT::kPrintLHS) return;
		
	fOut << "\nLHS matrix:\n\n";
	fOut << fMatrix << "\n\n";
}
