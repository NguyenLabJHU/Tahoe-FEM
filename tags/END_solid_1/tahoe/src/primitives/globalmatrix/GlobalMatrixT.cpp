/* $Id: GlobalMatrixT.cpp,v 1.5 2001-05-09 17:31:38 paklein Exp $ */
/* created: paklein (03/23/1997)                                          */
/* Virtual base class for all global matrix objects                       */

#include "GlobalMatrixT.h"
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include "Constants.h"
#include "dArrayT.h"

/* cconstructor */
GlobalMatrixT::GlobalMatrixT(ostream& out, int check_code):
	fOut(out),
	fCheckCode(check_code),
	fLocNumEQ(0),	
	fTotNumEQ(0),
	fStartEQ(0),
	fIsFactorized(0)
{
	if (fCheckCode < kNoCheck ||
	    fCheckCode > kPrintSolution) throw eBadInputValue;
}

GlobalMatrixT::~GlobalMatrixT(void) { }

/*
* Set the diagonal position matrix and allocate space.
* for the matrix.
*/
void GlobalMatrixT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	/* set dimensions */
	fTotNumEQ = tot_num_eq;
	fLocNumEQ = loc_num_eq;
	fStartEQ  = start_eq;

	/* output */
	fOut << "\n E q u a t i o n    S y s t e m    D a t a :\n\n";
	fOut << " Local number of equations . . . . . . . . . . . = " << fLocNumEQ << '\n';
	fOut << " Total number of equations . . . . . . . . . . . = " << fTotNumEQ << '\n';

	/* consistency */
	if (fLocNumEQ > fTotNumEQ)
	{
		cout << "\n GlobalMatrixT::Initialize: local number of equations " << fLocNumEQ << '\n'
		     <<   "     cannot be greater than the total number of equations "
		     << fTotNumEQ << endl;
		throw eGeneralFail;
	}

	/* must have at least 1 active equation */
	if (fLocNumEQ < 1)
	{
		cout << "\n GlobalMatrixT::Initialize: WARNING: no active equations"
		     << endl;
	}
	else if (fStartEQ < 1) /* active equation numbers must be > 0 */
	{
		cout << "\n GlobalMatrixT::Initialize: active equation must be > 0" << endl;
		throw eGeneralFail;
	}
}

/* set all matrix values to 0.0 */
void GlobalMatrixT::Clear(void) { fIsFactorized = 0; }

/*
* Solve the system for the vector given, returning the result
* in the same array
*/
void GlobalMatrixT::Solve(dArrayT& result)
{
	if (!fIsFactorized)
	{
		/* store original precision */
		int old_precision = fOut.precision();
	
		/* rank checks before factorization */
		fOut.precision(12);
		PrintLHS();
		fOut.precision(old_precision);
	
		/* factorize */
		fIsFactorized = 1;
		Factorize();
		
		/* rank checks after factorization */
		fOut.precision(12);
		PrintZeroPivots();
		PrintAllPivots();
		fOut.precision(old_precision);
	}

	/* output before solution */
	PrintRHS(result);

	/* find new search direction */
	BackSubstitute(result);

	/* output after solution */
	PrintSolution(result);
	
//TEMP: right result of linear solve	
//cout << "\n Writing update to output file: update.data" << endl;
//ofstream out("update.data");
//out.precision(12);
//out.setf(ios::showpoint);
//out.setf(ios::right, ios::adjustfield);
//out.setf(ios::scientific, ios::floatfield);
//out << result << endl;
}

/* strong manipulation functions 
 * NOTE: These must be overridden to provide support for these functions.
 *       By default, these all throw exceptions. These could be pure
 *       virtual, but that requires updating all derived matrix types */
void GlobalMatrixT::OverWrite(const ElementMatrixT& elMat, const nArrayT<int>& eqnos)
{
#pragma unused(elMat)
#pragma unused(eqnos)
	cout << "\n GlobalMatrixT::OverWrite: not implemented" << endl;
	throw eGeneralFail;
}

void GlobalMatrixT::Disassemble(dMatrixT& elMat, const nArrayT<int>& eqnos) const
{
#pragma unused(elMat)
#pragma unused(eqnos)
	cout << "\n GlobalMatrixT::Disassemble: not implemented" << endl;
	throw eGeneralFail;
}

void GlobalMatrixT::DisassembleDiagonal(dArrayT& diagonals, const nArrayT<int>& eqnos) const
{
#pragma unused(diagonals)
#pragma unused(eqnos)
	cout << "\n GlobalMatrixT::DisassembleDiagonal: not implemented" << endl;
	throw eGeneralFail;
}

/* assignment operator */
GlobalMatrixT& GlobalMatrixT::operator=(const GlobalMatrixT& RHS)
{
	fCheckCode    = RHS.fCheckCode;
	fLocNumEQ     = RHS.fLocNumEQ;
	fTotNumEQ     = RHS.fTotNumEQ;
	fStartEQ      = RHS.fStartEQ;
	fIsFactorized = RHS.fIsFactorized;
	return *this;
}

/**************************************************************************
* Protected
**************************************************************************/

void GlobalMatrixT::PrintRHS(const dArrayT& RHS) const
{
	if (fCheckCode != kPrintRHS) return;
	
	/* increase output stream precision */
	double* p = RHS.Pointer();

	int high_precision = 12;
	fOut.precision(high_precision);
	int d_width = OutputWidth(fOut, p);

	fOut << "\n RHS vector:\n\n";
	fOut << setw(kIntWidth)    << "eqn no.";
	fOut << setw(d_width) << "RHS\n\n";
	for (int i = 0; i < fLocNumEQ; i++)
	{
		fOut << setw(kIntWidth) << i + 1;
		fOut << setw(d_width) << *p++ << '\n';
	}
	
	fOut << endl;

	/* restore stream precision */
	fOut.precision(kPrecision);
}

void GlobalMatrixT::PrintSolution(const dArrayT& solution) const
{
	if (fCheckCode != kPrintSolution) return;
	
	/* increase output stream precision */
	double* p = solution.Pointer();

	int high_precision = 12;
	fOut.precision(high_precision);
	int d_width = OutputWidth(fOut, p);

	fOut << "\n solution vector:\n\n";
	fOut << setw(kIntWidth)    << "eqn no.";
	fOut << setw(d_width) << "solution\n\n";

	for (int i = 0; i < fLocNumEQ; i++)
	{
		fOut << setw(kIntWidth) << i + 1;
		fOut << setw(d_width) << *p++ << '\n';
	}
	
	fOut << endl;

	/* restore stream precision */
	fOut.precision(kPrecision);
}
