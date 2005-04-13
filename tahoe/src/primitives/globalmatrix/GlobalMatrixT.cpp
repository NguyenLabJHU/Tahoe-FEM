/* $Id: GlobalMatrixT.cpp,v 1.22 2005-04-13 17:40:37 paklein Exp $ */
/* created: paklein (03/23/1997) */
#include "GlobalMatrixT.h"
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include "toolboxConstants.h"
#include "dArrayT.h"

using namespace Tahoe;

/* initialize static data */
int GlobalMatrixT::sOutputCount = 0;

/* cconstructor */
GlobalMatrixT::GlobalMatrixT(ostream& out, int check_code):
	fOut(out),
	fCheckCode(check_code),
	fLocNumEQ(0),	
	fTotNumEQ(0),
	fStartEQ(0)
{

}

GlobalMatrixT::GlobalMatrixT(const GlobalMatrixT& source):
	fOut(source.fOut),
	fCheckCode(kNoCheck),
	fLocNumEQ(0),	
	fTotNumEQ(0),
	fStartEQ(0)
{
	GlobalMatrixT::operator=(source);
}

GlobalMatrixT::~GlobalMatrixT(void) { }

/*
* Set the diagonal position matrix and allocate space.
* for the matrix.
*/
void GlobalMatrixT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	const char caller[] = "GlobalMatrixT::Initialize";

	/* set dimensions */
	fTotNumEQ = tot_num_eq;
	fLocNumEQ = loc_num_eq;
	fStartEQ  = start_eq;

	/* consistency */
	if (fLocNumEQ > fTotNumEQ)
		ExceptionT::GeneralFail(caller, "local number of equations %d cannot exceed the total number of equations %d", 
			fLocNumEQ, fTotNumEQ);

	/* active equation numbers must be > 0 */
	if (fStartEQ < 1) ExceptionT::GeneralFail(caller, "active equation numbers must be > 0");
}

/* write information to output stream */
void GlobalMatrixT::Info(ostream& out) {
	out << "\n E q u a t i o n    S y s t e m    D a t a :\n\n";
	out << " Local number of equations . . . . . . . . . . . = " << fLocNumEQ << '\n';
	out << " Total number of equations . . . . . . . . . . . = " << fTotNumEQ << '\n';
	out.flush();
}

/*
* Solve the system for the vector given, returning the result
* in the same array
*/
bool GlobalMatrixT::Solve(dArrayT& result)
{
	/* catch any exceptions */
	try 
	{
		/* store original precision */
		int old_precision = fOut.precision();
		fOut.precision(12);
	
		/* rank checks before factorization */
		PrintLHS();
	
		/* factorize/precondition */
		Factorize();
	
		/* rank checks after factorization */
		PrintZeroPivots();
		PrintAllPivots();

		/* output before solution */
		PrintRHS(result);

		/* find new search direction */
		BackSubstitute(result);

		/* output after solution */
		PrintSolution(result);

		/* restore precision */
		fOut.precision(old_precision);
	}
	catch (ExceptionT::CodeT error) {
		cout << "\n GlobalMatrixT::Solve: caught exception: " << error << endl;
		return false;
	}
	return true;
}

/* strong manipulation functions 
 * NOTE: These must be overridden to provide support for these functions.
 *       By default, these all throw ExceptionT::xceptions. These could be pure
 *       virtual, but that requires updating all derived matrix types */
void GlobalMatrixT::OverWrite(const ElementMatrixT& elMat, const nArrayT<int>& eqnos)
{
#pragma unused(elMat)
#pragma unused(eqnos)
	ExceptionT::GeneralFail("GlobalMatrixT::OverWrite", "not implemented");
}

void GlobalMatrixT::Disassemble(dMatrixT& elMat, const nArrayT<int>& eqnos) const
{
#pragma unused(elMat)
#pragma unused(eqnos)
	ExceptionT::GeneralFail("GlobalMatrixT::Disassemble", "not implemented");
}

void GlobalMatrixT::DisassembleDiagonal(dArrayT& diagonals, const nArrayT<int>& eqnos) const
{
#pragma unused(diagonals)
#pragma unused(eqnos)
	ExceptionT::GeneralFail("GlobalMatrixT::DisassembleDiagonal", "not implemented");
}

/* compute the sum of the elements on the prescribed row/col */
double GlobalMatrixT::AbsRowSum(int rownum) const {
#pragma unused(rownum)
	ExceptionT::GeneralFail("GlobalMatrixT::AbsRowSum", "not implemented");
	return 0.0;
}

/* assignment operator */
GlobalMatrixT& GlobalMatrixT::operator=(const GlobalMatrixT& rhs)
{
	fCheckCode = rhs.fCheckCode;
	fLocNumEQ  = rhs.fLocNumEQ;
	fTotNumEQ  = rhs.fTotNumEQ;
	fStartEQ   = rhs.fStartEQ;

	return *this;
}

/* matrix-vector product */
void GlobalMatrixT::Multx(const dArrayT& x, dArrayT& b) const 
{ 
#pragma unused(x)
#pragma unused(b)
	ExceptionT::GeneralFail("GlobalMatrixT::Multx", "not implemented");
}

/* Tranpose[matrix]-vector product */
void GlobalMatrixT::MultTx(const dArrayT& x, dArrayT& b) const 
{
#pragma unused(x)
#pragma unused(b)
	ExceptionT::GeneralFail("GlobalMatrixT::MultTx", "not implemented");
}

/* vector-matrix-vector product */
double GlobalMatrixT::MultmBn(const dArrayT& m, const dArrayT& n) const
{
#pragma unused(m)
#pragma unused(n)
	ExceptionT::GeneralFail("GlobalMatrixT::MultmBn", "not implemented");
	return 0;
}

/* return the values along the diagonal of the matrix */
bool GlobalMatrixT::CopyDiagonal(dArrayT& diags) const
{
#pragma unused(diags)
	return false;
}

/**************************************************************************
 * Protected
 **************************************************************************/

void GlobalMatrixT::PrintRHS(const dArrayT& RHS) const
{
	if (fCheckCode != kPrintRHS) return;
	
	/* increase output stream precision */
	const double* p = RHS.Pointer();

	int high_precision = 12;
	fOut.precision(high_precision);
	int d_width = OutputWidth(fOut, p);

	fOut << "\n RHS vector:\n\n";
	fOut << setw(kIntWidth)    << "loc eq."
	     << setw(kIntWidth)    << "glb eq." 
	     << setw(d_width)      << "RHS\n\n";
	for (int i = 0; i < fLocNumEQ; i++)
		fOut << setw(kIntWidth) << i + 1
		     << setw(kIntWidth) << fStartEQ + i
		     << setw(d_width) << *p++ << '\n';	
	fOut << endl;

	/* restore stream precision */
	fOut.precision(kPrecision);
}

void GlobalMatrixT::PrintSolution(const dArrayT& solution) const
{
	if (fCheckCode != kPrintSolution) return;
	
	/* increase output stream precision */
	const double* p = solution.Pointer();

	int high_precision = 12;
	fOut.precision(high_precision);
	int d_width = OutputWidth(fOut, p);

	fOut << "\n solution vector:\n\n";
	fOut << setw(kIntWidth)    << "loc eq."
	     << setw(kIntWidth)    << "glb eq." 
	     << setw(d_width)      << "solution\n\n";
	for (int i = 0; i < fLocNumEQ; i++)
		fOut << setw(kIntWidth) << i + 1
		     << setw(kIntWidth) << fStartEQ + i
		     << setw(d_width) << *p++ << '\n';	
	fOut << endl;

	/* restore stream precision */
	fOut.precision(kPrecision);
}
