/* $Id: SPOOLESMatrixT_mpi.cpp,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (09/13/2000)                                          */
/* SPOOLES matrix solver                                                  */

#include "SPOOLESMatrixT_mpi.h"

/* library support options */
#ifdef __SPOOLES__
#ifdef __MPI__
#include "SPOOLESMPI.h"
#include "StringT.h"
#include "MSRBuilderT.h"

/* message file name */
const char SPOOLES_FILE_ROOT[] = "SPOOLES";
const char  SPOOLES_FILE_EXT[] = ".out";

/* constuctor */
SPOOLESMatrixT_mpi::SPOOLESMatrixT_mpi(ostream& out, int check_code,
	bool symmetric, bool pivoting):
	SPOOLESMatrixT(out, check_code, symmetric, pivoting)
{

}

/* assignment operator */
GlobalMatrixT& SPOOLESMatrixT_mpi::operator=(const GlobalMatrixT& RHS)
{
#pragma unused(RHS)

	cout << "\n SPOOLESMatrixT_mpi::operator=: not implemented" << endl;
	throw eGeneralFail;
}

/*************************************************************************
* Protected
*************************************************************************/

/* determine new search direction and put the results in result */
void SPOOLESMatrixT_mpi::BackSubstitute(dArrayT& result)
{
	/* flag should not be set */
	if (fIsFactorized) throw eGeneralFail;

#if 0
ofstream out("rcv.out");
out << "graph:\n";
fMSRBuilder->Write(out);
out << "\nMSR data:\n";
fMSRBuilder->WriteMSRData(out, fupdate, fbindx);
#endif

	/* convert matrix to RCV */
	iArrayT r, c;
	dArrayT v;
	GenerateRCV(r, c, v);

//DEBUG
#if 0
//ofstream out("rcv.out");
for (int i = 0; i < r.Length(); i++)
	out << setw(kIntWidth) << r[i]
	    << setw(kIntWidth) << c[i]
		<< setw(kDoubleWidth) << v[i] << endl;
cout << "\n SPOOLESMatrixT::BackSubstitute: EXIT" << endl;
throw;
#endif
//DEBUG

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
			throw eGeneralFail;
		}
	}
	else
	{
		MPI_Comm comm = MPI_COMM_WORLD; //TEMP

		/* process rank */
		int rank;
		if (MPI_Comm_rank(comm, &rank) != MPI_SUCCESS) throw eMPIFail;

		/* message file name */
		StringT spooles_file(SPOOLES_FILE_ROOT);
		spooles_file.Append(".p", rank);
		spooles_file.Append(SPOOLES_FILE_EXT);

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
			throw eGeneralFail;
		}
	}
}
#endif /* __MPI__ */
#endif /* __SPOOLES__ */
