/* $Id: SPOOLESMatrixT_mpi.cpp,v 1.14 2004-03-14 07:35:40 paklein Exp $ */
/* created: paklein (09/13/2000) */

#include "SPOOLESMatrixT_mpi.h"

/* library support options */
#ifdef __SPOOLES_MPI__
#ifdef __TAHOE_MPI__
#include "SPOOLESMPI.h"
#include "StringT.h"
#include "MSRBuilderT.h"
#include "CommunicatorT.h"

using namespace Tahoe;

/* message file name */
const char SPOOLES_FILE_ROOT[] = "SPOOLES";
const char  SPOOLES_FILE_EXT[] = ".out";

/* constuctor */
SPOOLESMatrixT_mpi::SPOOLESMatrixT_mpi(ostream& out, int check_code,
	bool symmetric, bool pivoting, CommunicatorT& comm):
	SPOOLESMatrixT(out, check_code, symmetric, pivoting),
	fComm(comm)
{

}

/*************************************************************************
* Protected
*************************************************************************/

/* determine new search direction and put the results in result */
void SPOOLESMatrixT_mpi::BackSubstitute(dArrayT& result)
{
	const char caller[] = "SPOOLESMatrixT_mpi::BackSubstitute";

	/* convert matrix to RCV */
	iArrayT r, c;
	dArrayT v;
	GenerateRCV(r, c, v, 1.0e-15);

	/* serial driver provided by in SPOOLES documentation */
	int msglvl = 0; //  0: nothing
	                //  1: scalar output (timing data) only
	                // >1: verbose
	int matrix_type = SPOOLES_REAL;
	int symmetry_flag = (fSymmetric) ? SPOOLES_SYMMETRIC : SPOOLES_NONSYMMETRIC;
	int pivoting_flag = (fPivoting) ? SPOOLES_PIVOTING : SPOOLES_NO_PIVOTING;
	int seed = 1;

	if (fLocNumEQ == fTotNumEQ)
	{
		StringT spooles_file(SPOOLES_FILE_ROOT);
		spooles_file.Append(SPOOLES_FILE_EXT);
		int OK = LU_serial_driver(msglvl, spooles_file, matrix_type, symmetry_flag,
			pivoting_flag, seed, result.Length(), result.Pointer(),
			r.Length(), r.Pointer(), c.Pointer(), v.Pointer());
		if (OK != 1) 
			ExceptionT::BadJacobianDet(caller, "LU_serial_driver returned: %d", OK);
	}
	else
	{
		/* message file name */
		StringT spooles_file(SPOOLES_FILE_ROOT);
		spooles_file.Append(".p", fComm.Rank());
		spooles_file.Append(SPOOLES_FILE_EXT);
		
		//TEMP - SPOOLES v2.2 does not seem to solve non-symmetric matricies with pivoting disabled
		if (!fSymmetric && !fPivoting)
			cout << "\n SPOOLESMatrixT_mpi::BackSubstitute: WARNING: SPOOLES v2.2 does not solve\n"
			     <<   "     non-symmetric systems correctly with pivoting disabled."<< endl;

		/* driver */
#ifndef __MWERKS__
		
		/* the MPI_Comm */
		MPI_Comm comm = fComm;
		int OK = LU_MPI_driver(msglvl, spooles_file, matrix_type, symmetry_flag,
			pivoting_flag, seed, fTotNumEQ, result.Length(), fupdate.Pointer(),
			result.Pointer(), r.Length(), r.Pointer(), c.Pointer(), v.Pointer(), &comm);
#else
		int OK = 0;
		cout << "\n SPOOLESMatrixT_mpi::BackSubstitute: MPI SPOOLES not supported by MacMPI" << endl;
#endif /* __MWERKS__ */
		if (OK != 1)
			ExceptionT::BadJacobianDet(caller, "LU_MPI_driver returned: %d", OK);
	}
}
#endif /* __TAHOE_MPI__ */
#endif /* __SPOOLES_MPI__ */
