/* $Id: TrilinosAztecT.cpp,v 1.1 2006-11-01 05:15:31 paklein Exp $ */
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

	/* set-up Epetra structures */
#ifdef __TAHOE_MPI__
	Epetra_MpiComm epetra_comm(fComm);
#else
	Epetra_SerialComm epetra_comm;
#endif
	iArrayT active_tmp;
	active_tmp.Alias(factive);
	active_tmp--; /* solver uses 0 indexing */
	Epetra_Map map(fTotNumEQ, fLocNumEQ, active_tmp.Pointer(), 0, epetra_comm);
	iArrayT row_count(fLocNumEQ);
	for (int i = 0; i < fLocNumEQ - 1; i++) {
		row_count[i] = frowptr[i+1] - frowptr[i];
	}
	row_count[fLocNumEQ-1] = fnzval.Length() - frowptr[fLocNumEQ-1];
	Epetra_CrsMatrix A(Copy, map, row_count.Pointer(), true);

	/* copy data into epetra_matrix */
	for (int i = 0; i < factive.Length(); i++) {
		int offset = frowptr[i];
		int ret = A.InsertGlobalValues(factive[i], row_count[i], fnzval.Pointer(offset), fcolind.Pointer(offset));
		if (ret != 0) {
			ExceptionT::GeneralFail(caller, "InsertGlobalValues: error %d in row %d", ret, i+1);
		}
	}  
	int ret = A.FillComplete();
	if (ret != 0) {
		ExceptionT::GeneralFail(caller, "FillComplete: error %d", ret);
	}

	/* Create x and b vectors */
	Epetra_Vector b(Copy, map, result.Pointer());
	Epetra_Vector x(View, map, result.Pointer());

	/* create linear problem */
	Epetra_LinearProblem problem(&A, &x, &b);

	/* create AztecOO instance */
	AztecOO solver(problem);

	/* set options */
	solver.SetAztecOption(AZ_precond, AZ_Jacobi);

	/* solve */
	solver.Iterate(1000, 1.0E-8);
	active_tmp++; /* restore first active row = 1 */
}

#endif /* __TRILINOS__ */
