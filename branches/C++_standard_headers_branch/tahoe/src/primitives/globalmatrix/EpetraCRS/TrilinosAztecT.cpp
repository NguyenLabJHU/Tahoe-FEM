/* $Id: TrilinosAztecT.cpp,v 1.2 2006-11-14 04:31:40 paklein Exp $ */
#include "TrilinosAztecT.h"

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
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#ifdef __TAHOE_MPI__
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "AztecOO.h"

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

using namespace Tahoe;

/* constructor */
TrilinosAztecT::TrilinosAztecT(ostream& out, int check_code, const CommunicatorT& comm):
	EpetraCRSMatrixT(out, check_code, comm)
{

}

/* copy constructor */
TrilinosAztecT::TrilinosAztecT(const TrilinosAztecT& rhs):
	EpetraCRSMatrixT(rhs)
{
	TrilinosAztecT::operator=(rhs);
}

TrilinosAztecT& TrilinosAztecT::operator=(const TrilinosAztecT&)
{
	ExceptionT::GeneralFail("TrilinosAztecT::operator=", "not implemented");
	return *this;
}

/* return a clone of self */
GlobalMatrixT* TrilinosAztecT::Clone(void) const {
	return new TrilinosAztecT(*this);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* solution driver */
void TrilinosAztecT::BackSubstitute(dArrayT& result)
{
	const char caller[] = "TrilinosAztecT::BackSubstitute";

	/* get Epetra matrix */
	Epetra_CrsMatrix* A = Translate();

	/* Create x and b vectors */
	Epetra_Vector b(Copy, *fepetra_map, result.Pointer());
	Epetra_Vector x(View, *fepetra_map, result.Pointer());

	/* create linear problem */
	Epetra_LinearProblem problem(A, &x, &b);

	/* create AztecOO instance */
	AztecOO solver(problem);

	/* set options */
	solver.SetAztecOption(AZ_precond, AZ_Jacobi);

	/* solve */
	solver.Iterate(1000, 1.0E-8);
	
	/* clean-up */
	delete A;
}

#endif /* __TRILINOS__ */
