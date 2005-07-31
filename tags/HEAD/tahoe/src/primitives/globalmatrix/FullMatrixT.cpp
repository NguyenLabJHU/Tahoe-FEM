/* $Id: FullMatrixT.cpp,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (03/07/1998)                                          */
/* Virtual base class for all global matrix objects                       */

#include "FullMatrixT.h"
#include <iostream.h>
#include <iomanip.h>
#include "Constants.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "ElementMatrixT.h"

/* constructor */
FullMatrixT::FullMatrixT(ostream& out,int check_code):
	GlobalMatrixT(out, check_code)
{

}

/* set the internal matrix structure.
* NOTE: do not call Initialize() equation topology has been set
* with AddEquationSet() for all equation sets */
void FullMatrixT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	/* inherited */
	GlobalMatrixT::Initialize(tot_num_eq, loc_num_eq, start_eq);

	/* check */
	if (tot_num_eq != loc_num_eq)
	{
		cout << "\n FullMatrixT::Initialize: expecting total number of equations\n"
		     <<   "     " << tot_num_eq
		     << " to be equal to the local number of equations " << loc_num_eq << endl;
		throw eGeneralFail;
	}
	
	/* allocate work space */
	fMatrix.Allocate(fLocNumEQ);
}

/* set all matrix values to 0.0 */
void FullMatrixT::Clear(void)
{
	/* inherited */
	GlobalMatrixT::Clear();
	
	/* clear values */
	fMatrix = 0.0;	
}

/* add element group equations to the overall topology.
* NOTE: assembly positions (equation numbers) = 1...fNumEQ
* equations can be of fixed size (iArray2DT) or
* variable length (RaggedArray2DT) */
void FullMatrixT::AddEquationSet(const iArray2DT& eqset)
{
#pragma unused(eqset)
// no equation data needed
}

void FullMatrixT::AddEquationSet(const RaggedArray2DT<int>& eqset)
{
#pragma unused(eqset)
// no equation data needed
}

/* assemble the element contribution into the LHS matrix - assumes
* that elMat is square (n x n) and that eqnos is also length n.
* NOTE: assembly positions (equation numbers) = 1...fNumEQ */
void FullMatrixT::Assemble(const ElementMatrixT& elMat, const iArrayT& eqnos)
{
#if __option(extended_errorcheck)
	/* dimension checking */
	if (elMat.Rows() != eqnos.Length() ||
	    elMat.Cols() != eqnos.Length()) throw eSizeMismatch;
#endif

	/* element matrix format */
	ElementMatrixT::FormatT format = elMat.Format();

	if (format == ElementMatrixT::kDiagonal)
	{
		/* from diagonal only */
		double* pelMat = elMat.Pointer();
		int inc = elMat.Rows() + 1;
	
		int*    peq    = eqnos.Pointer();		
		for (int i = 0; i < elMat.Length(); i++)
		{
			int eq = *peq++;
			
			/* active dof's only */
			if (eq-- > 0) fMatrix(eq,eq) += *pelMat;
			
			pelMat += inc;
		}
	}
	else
	{
		/* copy to full symmetric */
		if (elMat.Format() == ElementMatrixT::kSymmetricUpper) elMat.CopySymmetric();

		int*    peq    = eqnos.Pointer();	
		for (int col = 0; col < elMat.Cols(); col++)
		{
			int        eqc = eqnos[col];
			double* pelMat = elMat(col);
		
			/* active dof's only */
			if (eqc-- > 0)
			{		
				int* peqr = eqnos.Pointer();
						
				for (int row = 0; row < elMat.Rows(); row++)
				{
					int eqr = *peqr++;
				
					/* active dof's only */
					if (eqr-- > 0) fMatrix(eqr,eqc) += *pelMat;
					
					pelMat++;
				}
			}
		}
	}
}

/* strong manipulation functions */
void FullMatrixT::OverWrite(const ElementMatrixT& elMat, const iArrayT& eqnos)
{
#if __option(extended_errorcheck)
	/* dimension checking */
	if (elMat.Rows() != eqnos.Length() ||
	    elMat.Cols() != eqnos.Length()) throw eSizeMismatch;
#endif

	/* copy to full symmetric */
	if (elMat.Format() == ElementMatrixT::kSymmetricUpper) elMat.CopySymmetric();

	int* peq = eqnos.Pointer();	
	for (int col = 0; col < elMat.Cols(); col++)
	{
		int        eqc = eqnos[col];
		double* pelMat = elMat(col);
		
		/* active dof's only */
		if (eqc-- > 0)
		{		
			int* peqr = eqnos.Pointer();
					
			for (int row = 0; row < elMat.Rows(); row++)
			{
				int eqr = *peqr++;
			
				/* active dof's only */
				if (eqr-- > 0) fMatrix(eqr,eqc) += *pelMat;
				
				pelMat++;
			}
		}
	}
}

void FullMatrixT::Disassemble(dMatrixT& elMat, const iArrayT& eqnos) const
{
#if __option(extended_errorcheck)
	/* dimension checking */
	if (elMat.Rows() != eqnos.Length() ||
	    elMat.Cols() != eqnos.Length()) throw eSizeMismatch;
#endif

	int* peq = eqnos.Pointer();
	
	for (int col = 0; col < elMat.Cols(); col++)
	{
		int        eqc = eqnos[col];
		double* pelMat = elMat(0);
		
		/* active dof's only */
		if (eqc-- > 0)
		{		
			int* peqr = eqnos.Pointer();
					
			for (int row = 0; row < elMat.Rows(); row++)
			{
				int eqr = *peqr++;
			
				/* active dof's only */
				if (eqr-- > 0)
					*pelMat = fMatrix(eqr,eqc);
				else
					*pelMat = 0.0;
				
				pelMat++;
			}
		}
		else
			/* clear column */
			elMat.SetCol(col, 0.0);
	}
}

/* assignment operator */
GlobalMatrixT& FullMatrixT::operator=(const GlobalMatrixT& RHS)
{
	/* inherited */
	GlobalMatrixT::operator=(RHS);

#ifdef __NO_RTTI__
	const FullMatrixT* pRHS = (const FullMatrixT*) (&RHS);
#else
	const FullMatrixT* pRHS = dynamic_cast<const FullMatrixT*>(&RHS);
	if (!pRHS) throw eGeneralFail;
#endif

	/* copy matrix */
	fMatrix = pRHS->fMatrix;   	

	return *this;
}

/* number scope and reordering */
GlobalMatrixT::EquationNumberScopeT FullMatrixT::EquationNumberScope(void) const
{
	return kLocal;
}

bool FullMatrixT::RenumberEquations(void) const { return false; }

/**************************************************************************
* Protected
**************************************************************************/

/* precondition matrix */
void FullMatrixT::Factorize(void)
{
//TEMP: no full, nonsymmetric factorization implemented
	fIsFactorized = 0; // undo GlobalMatrixT flag set
}
	
/* determine new search direction and put the results in result */
void FullMatrixT::BackSubstitute(dArrayT& result)
{
//TEMP: no full, nonsymmetric factorization implemented
	if (fIsFactorized) throw eGeneralFail;

	fMatrix.LinearSolve(result);
}

/* rank check functions */
void FullMatrixT::PrintAllPivots(void) const
{
//TEMP: no full, nonsymmetric factorization implemented
}

void FullMatrixT::PrintZeroPivots(void) const
{
//TEMP: no full, nonsymmetric factorization implemented
}

void FullMatrixT::PrintLHS(void) const
{
	if (fCheckCode != GlobalMatrixT::kPrintLHS) return;
		
	fOut << "\nLHS matrix:\n\n";
	fOut << fMatrix << "\n\n";
}
