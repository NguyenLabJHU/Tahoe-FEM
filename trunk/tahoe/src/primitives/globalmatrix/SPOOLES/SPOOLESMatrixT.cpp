/* $Id: SPOOLESMatrixT.cpp,v 1.24 2005-01-07 21:23:02 paklein Exp $ */
/* created: paklein (09/13/2000) */
#include "SPOOLESMatrixT.h"

/* library support options */
#ifdef __SPOOLES__

#include "ElementMatrixT.h"
#include "LU_serial_driver.h"

#ifdef __SPOOLES_MT__
#include "SPOOLESMT.h"
#endif

using namespace Tahoe;

/* log file */
const char SPOOLES_FILE[] = "SPOOLES.out";

/* constuctor */
SPOOLESMatrixT::SPOOLESMatrixT(ostream& out, int check_code,
	bool symmetric, bool pivoting):
	MSRMatrixT(out, check_code, symmetric),
	fPivoting(pivoting),
	pLU_dat(NULL),
	fIsFactorized(false)
{

}

SPOOLESMatrixT::SPOOLESMatrixT(const SPOOLESMatrixT& source):
	MSRMatrixT(source),
	pLU_dat(NULL),
	fIsFactorized(false)	
{
	ExceptionT::GeneralFail("SPOOLESMatrixT::SPOOLESMatrixT", "not implemented");
}

/* destructor */
SPOOLESMatrixT::~SPOOLESMatrixT(void)
{
	if (pLU_dat && !LU_serial_driver_free(&pLU_dat))
		ExceptionT::GeneralFail("SPOOLESMatrixT::~SPOOLESMatrixT",
			"error freeing SPOOLES data");
	pLU_dat = NULL;
}

/* clear values for next assembly */
void SPOOLESMatrixT::Clear(void)
{
	const char caller[] = "SPOOLESMatrixT::Clear";

	/* inherited */
	MSRMatrixT::Clear();

	/* check - must be serial */
	if (fTotNumEQ != fLocNumEQ)
		ExceptionT::GeneralFail(caller, "total equations (%d) != local equations (%d)",
			fTotNumEQ, fLocNumEQ);

#ifndef __SPOOLES_MT__
	/* delete existing data */
	if (pLU_dat && !LU_serial_driver_free(&pLU_dat))
		ExceptionT::GeneralFail(caller, "error freeing SPOOLES data");

	/* initialize new data */
	int matrix_type = SPOOLES_REAL;
	int symmetry_flag = (fSymmetric) ? SPOOLES_SYMMETRIC : SPOOLES_NONSYMMETRIC;
	int pivoting_flag = (fPivoting) ? SPOOLES_PIVOTING : SPOOLES_NO_PIVOTING;
	int seed = 1; /* any seed will do here */
	if (!LU_serial_driver_init(matrix_type, symmetry_flag, pivoting_flag, seed, fTotNumEQ, &pLU_dat))
		ExceptionT::GeneralFail(caller, "error initializing SPOOLES data");
#endif /* __SPOOLES_MT__ */

	/* reset flag */
	fIsFactorized = false;
}

/* assignment operator */
SPOOLESMatrixT& SPOOLESMatrixT::operator=(const SPOOLESMatrixT&)
{
	const char caller[] = "SPOOLESMatrixT::operator=";
	ExceptionT::GeneralFail(caller, "not implemented");
	return *this;
}

/** return a clone of self */
GlobalMatrixT* SPOOLESMatrixT::Clone(void) const {
	return new SPOOLESMatrixT(*this);
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* precondition matrix */
void SPOOLESMatrixT::Factorize(void)
{
	/* quick exit */
	if (fIsFactorized) return;

#ifndef __SPOOLES_MT__
	/* convert matrix to RCV */
	iArrayT r, c;
	dArrayT v;
	GenerateRCV(r, c, v, 1.0e-15);

	/* logging level:
		 0: nothing
		 1: scalar output (timing data) only
		>1: verbose */
	int msglvl = 0; 

	/* compute factorization */
	int OK = LU_serial_driver_factorize(msglvl, SPOOLES_FILE,
		r.Length(), r.Pointer(), c.Pointer(), v.Pointer(), 
		&pLU_dat);

	if (OK != 1) ExceptionT::BadJacobianDet("SPOOLESMatrixT::Factorize", "LU driver returned %d", OK);
#endif /* __SPOOLES_MT__ */

	/* set flag */
	fIsFactorized = true;
}

/* determine new search direction and put the results in result */
void SPOOLESMatrixT::BackSubstitute(dArrayT& result)
{
	const char caller[] = "SPOOLESMatrixT::BackSubstitute";
 
 	/* check */
	if (!fIsFactorized) ExceptionT::GeneralFail(caller, "matrix is not factorized");

	 int msglvl = 0; //  0: nothing
	 //  1: scalar output (timing data) only
	 // >1: verbose
	int OK;
 
#ifdef __SPOOLES_MT__
	/* convert matrix to RCV */
	iArrayT r, c;
	dArrayT v;
	GenerateRCV(r, c, v, 1.0e-15);

	int matrix_type = SPOOLES_REAL;
	int symmetry_flag = (fSymmetric) ? SPOOLES_SYMMETRIC : SPOOLES_NONSYMMETRIC;
	int pivoting_flag = (fPivoting) ? SPOOLES_PIVOTING : SPOOLES_NO_PIVOTING;
	int seed = 1; /* any seed will do here */

	/* call MT driver, __NUM_THREADS__ comes from make macro */
	OK = LU_MT_driver(msglvl, SPOOLES_FILE, matrix_type, symmetry_flag,
		pivoting_flag, seed, result.Length(), result.Pointer(),
		r.Length(), r.Pointer(), c.Pointer(), v.Pointer(), __NUM_THREADS__);
#else /* not __SPOOLES_MT__ */

	/* back substitute */
	OK = LU_serial_driver_solve(msglvl, SPOOLES_FILE, result.Pointer(), &pLU_dat);

#endif /* __SPOOLES_MT__ */
 
	if (OK != 1) ExceptionT::BadJacobianDet(caller, "LU driver returned %d", OK);
}

#endif /* __SPOOLES__ */
