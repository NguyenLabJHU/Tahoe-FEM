/* $Id: SPOOLESMatrixT_mpi.cpp,v 1.20 2005-04-13 22:27:07 paklein Exp $ */
/* created: paklein (09/13/2000) */
#include "SPOOLESMatrixT_mpi.h"

/* library support options */
#ifdef __SPOOLES_MPI__
#ifdef __TAHOE_MPI__
#include "LU_MPI_driver.h"
//TEMP
#include "SPOOLESMPI.h"

#include "StringT.h"
#include "MSRBuilderT.h"
#include "CommunicatorT.h"

using namespace Tahoe;

/* message file name */
const char SPOOLES_FILE_ROOT[] = "SPOOLES";
const char  SPOOLES_FILE_EXT[] = ".out";

//old code or new code
#undef OLD_CODE

/* constuctor */
SPOOLESMatrixT_mpi::SPOOLESMatrixT_mpi(ostream& out, int check_code,
	bool symmetric, bool pivoting, int message_level, const CommunicatorT& comm):
	SPOOLESMatrixT(out, check_code, symmetric, pivoting, message_level, comm)
{

}

/* destructor */
SPOOLESMatrixT_mpi::~SPOOLESMatrixT_mpi(void)
{
	if (pLU_dat && !LU_MPI_driver_free(&pLU_dat))
		ExceptionT::GeneralFail("SPOOLESMatrixT_mpi::~SPOOLESMatrixT_mpi",
			"error freeing SPOOLES data");
	pLU_dat = NULL;
}

/* clear values for next assembly */
void SPOOLESMatrixT_mpi::Clear(void)
{
	const char caller[] = "SPOOLESMatrixT::Clear";

	/* inherited - don't call SPOOLESMatrixT::Clear */
	MSRMatrixT::Clear();

	/* system does not appear to be distributed */
	if (fLocNumEQ == fTotNumEQ)
		ExceptionT::GeneralFail(caller, "local and total number of equations the same: %d", fLocNumEQ);

	/* delete existing data */
	if (pLU_dat && !LU_MPI_driver_free(&pLU_dat))
		ExceptionT::GeneralFail(caller, "error freeing SPOOLES data");

	/* initialize new data */
	int matrix_type = SPOOLES_REAL;
	int symmetry_flag = (fSymmetric) ? SPOOLES_SYMMETRIC : SPOOLES_NONSYMMETRIC;
	int pivoting_flag = (fPivoting) ? SPOOLES_PIVOTING : SPOOLES_NO_PIVOTING;
	int seed = 1; /* anything will do */
	if (!LU_MPI_driver_init(matrix_type, symmetry_flag, pivoting_flag, seed, fTotNumEQ, fLocNumEQ, &pLU_dat))
		ExceptionT::GeneralFail(caller, "error initializing SPOOLES data");

	/* reset flag */
	fIsFactorized = false;
}

/* assignment operator */
SPOOLESMatrixT_mpi& SPOOLESMatrixT_mpi::operator=(const SPOOLESMatrixT_mpi&)
{
	ExceptionT::GeneralFail("SPOOLESMatrixT_mpi::operator=", "not implemented");
	return *this;
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* precondition matrix */
void SPOOLESMatrixT_mpi::Factorize(void)
{
#ifndef OLD_CODE
	/* quick exit */
	if (fIsFactorized) return;

	/* convert matrix to RCV */
	iArrayT r, c;
	dArrayT v;
	GenerateRCV(r, c, v, 1.0e-15);

	/* logging level:
		 0: nothing
		 1: scalar output (timing data) only
		>1: verbose */
	int msglvl = (fMessageLevel < 0) ? 0 : fMessageLevel; 

	/* message file name */
	StringT spooles_file(SPOOLES_FILE_ROOT);
	spooles_file.Append(".p", fComm.Rank());
	spooles_file.Append(SPOOLES_FILE_EXT);

	//TEMP - SPOOLES v2.2 does not seem to solve non-symmetric matricies with pivoting disabled
	if (!fSymmetric && !fPivoting)
		cout << "\n SPOOLESMatrixT_mpi::Factorize: WARNING: SPOOLES v2.2 does not solve\n"
		     <<   "     non-symmetric systems correctly with pivoting disabled."<< endl;

	/* the MPI_Comm */
	MPI_Comm comm = fComm;

	/* compute factorization */
	int OK = LU_MPI_driver_factorize(msglvl, spooles_file,
		r.Length(), r.Pointer(), c.Pointer(), v.Pointer(), 
		&pLU_dat, &comm);

	if (OK != 1) ExceptionT::BadJacobianDet("SPOOLESMatrixT_mpi::Factorize", "LU driver returned %d", OK);

	fIsFactorized = true;
#endif
}

/* determine new search direction and put the results in result */
void SPOOLESMatrixT_mpi::BackSubstitute(dArrayT& result)
{
	const char caller[] = "SPOOLESMatrixT_mpi::BackSubstitute";

#ifdef OLD_CODE
	/* convert matrix to RCV */
	iArrayT r, c;
	dArrayT v;
	GenerateRCV(r, c, v, 1.0e-15);

	/* serial driver provided by in SPOOLES documentation */
	int msglvl = (fMessageLevel < 0) ? 0 : fMessageLevel; 
	//  0: nothing
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

#else /* not OLD_CODE */

 	/* check */
	if (!fIsFactorized) ExceptionT::GeneralFail(caller, "matrix is not factorized");

	/* logging level:
		 0: nothing
		 1: scalar output (timing data) only
		>1: verbose */
	int msglvl = (fMessageLevel < 0) ? 0 : fMessageLevel;

	/* message file name */
	StringT spooles_file(SPOOLES_FILE_ROOT);
	spooles_file.Append(".p", fComm.Rank());
	spooles_file.Append(SPOOLES_FILE_EXT);
		
	/* the MPI_Comm */
	MPI_Comm comm = fComm;

	/* driver */
	int OK = LU_MPI_driver_solve(msglvl, spooles_file, fupdate.Pointer(), result.Pointer(), 
		&pLU_dat, &comm);

	/* check */
	if (OK != 1) ExceptionT::BadJacobianDet(caller, "LU_MPI_driver returned: %d", OK);
#endif /* OLD_CODE */
}
#endif /* __TAHOE_MPI__ */
#endif /* __SPOOLES_MPI__ */
