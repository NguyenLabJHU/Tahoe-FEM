/* $Id: SPOOLESMatrixT_mpi.cpp,v 1.9 2002-11-28 17:06:32 paklein Exp $ */
/* created: paklein (09/13/2000) */

#include "SPOOLESMatrixT_mpi.h"

/* library support options */
#ifdef __SPOOLES_MPI__
#ifdef __TAHOE_MPI__
#include "SPOOLESMPI.h"
#include "StringT.h"
#include "MSRBuilderT.h"

/* message file name */

using namespace Tahoe;

const char SPOOLES_FILE_ROOT[] = "SPOOLES";
const char  SPOOLES_FILE_EXT[] = ".out";

/* constuctor */
SPOOLESMatrixT_mpi::SPOOLESMatrixT_mpi(ostream& out, int check_code,
	bool symmetric, bool pivoting):
	SPOOLESMatrixT(out, check_code, symmetric, pivoting)
{

}

/*************************************************************************
* Protected
*************************************************************************/

/* determine new search direction and put the results in result */
void SPOOLESMatrixT_mpi::BackSubstitute(dArrayT& result)
{
	/* flag should not be set */
	if (fIsFactorized) throw ExceptionT::kGeneralFail;

	/* convert matrix to RCV */
	iArrayT r, c;
	dArrayT v;
	GenerateRCV(r, c, v, 1.0e-15);
	
	/* write matrix */
	if (fCheckCode == kPrintLHS) {
		int old_precision = fOut.precision();
		fOut.precision(12);
		int d_width = fOut.precision() + kDoubleExtra;
		fOut << "\n LHS matrix in {r,c,v} format:\n";
		fOut << " Number of values = " << r.Length() << '\n';
		for (int i = 0; i < r.Length(); i++)
			fOut << setw(kIntWidth) << r[i] + 1
			     << setw(kIntWidth) << c[i] + 1
			     << setw(d_width) << v[i] << '\n';
		fOut.flush();
		fOut.precision(old_precision);
	}

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
		{
			cout << "\n SPOOLESMatrixT_mpi::BackSubstitute: LU_serial_driver returned: "
			     << OK << endl;
			throw ExceptionT::kGeneralFail;
		}
	}
	else
	{
		MPI_Comm comm = MPI_COMM_WORLD; //TEMP

		/* process rank */
		int rank;
		if (MPI_Comm_rank(comm, &rank) != MPI_SUCCESS) throw ExceptionT::kMPIFail;

		/* message file name */
		StringT spooles_file(SPOOLES_FILE_ROOT);
		spooles_file.Append(".p", rank);
		spooles_file.Append(SPOOLES_FILE_EXT);
		
		//TEMP - SPOOLES v2.2 does not seem to solve non-symmetric matricies with pivoting disabled
		if (!fSymmetric && !fPivoting)
			cout << "\n SPOOLESMatrixT_mpi::BackSubstitute: WARNING: SPOOLES v2.2 does not solve\n"
			     <<   "     non-symmetric systems correctly with pivoting disabled."<< endl;

		/* driver */
#ifndef __MWERKS__
		int OK = LU_MPI_driver(msglvl, spooles_file, matrix_type, symmetry_flag,
			pivoting_flag, seed, fTotNumEQ, result.Length(), fupdate.Pointer(),
			result.Pointer(), r.Length(), r.Pointer(), c.Pointer(), v.Pointer(), &comm);
#else
		int OK = 0;
		cout << "\n SPOOLESMatrixT_mpi::BackSubstitute: MPI SPOOLES not supported by MacMPI" << endl;
#endif /* __MWERKS__ */
		if (OK != 1)
		{
			cout << "\n SPOOLESMatrixT_mpi::BackSubstitute: LU_MPI_driver returned: "
			     << OK << endl;
			throw ExceptionT::kGeneralFail;
		}
	}
}
#endif /* __TAHOE_MPI__ */
#endif /* __SPOOLES_MPI__ */
