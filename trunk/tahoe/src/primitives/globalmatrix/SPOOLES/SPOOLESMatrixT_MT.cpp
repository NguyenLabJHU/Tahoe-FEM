/* $Id: SPOOLESMatrixT_MT.cpp,v 1.4 2005-04-18 05:47:55 paklein Exp $ */
/* created: paklein (09/13/2000) */
#include "SPOOLESMatrixT_MT.h"

/* library support options */
#ifdef __SPOOLES_MT__

#include "ElementMatrixT.h"
#include "LU_MT_driver.h"

using namespace Tahoe;

/* log file */
const char SPOOLES_FILE[] = "SPOOLES.out";

//old code or new code
#undef OLD_CODE

/* constuctor */
SPOOLESMatrixT_MT::SPOOLESMatrixT_MT(ostream& out, int check_code,
	bool symmetric, bool pivoting, int message_level, int num_threads, const CommunicatorT& comm):
	MSRMatrixT(out, check_code, symmetric, comm),
	fPivoting(pivoting),
	fMessageLevel(message_level),
	pLU_dat(NULL),
	fIsFactorized(false),
	fNumThreads(num_threads)
{
/* must define thread model */
#if ! defined(__POSIX_THREADS__) && ! defined(__SOLARIS_THREADS__)
#error "requires __POSIX_THREADS__ or __SOLARIS_THREADS__"
ExceptionT::GeneralFail("SPOOLESMatrixT_MT::SPOOLESMatrixT_MT",
	"requires __POSIX_THREADS__ or __SOLARIS_THREADS__");
#endif

	/* check */
	if (fNumThreads < 2) 
		ExceptionT::GeneralFail("SPOOLESMatrixT_MT::SPOOLESMatrixT_MT",
			"expecting at least 2 threads %d", fNumThreads);
}

SPOOLESMatrixT_MT::SPOOLESMatrixT_MT(const SPOOLESMatrixT_MT& source):
	MSRMatrixT(source),
	fPivoting(source.fPivoting),
	fMessageLevel(source.fMessageLevel),
	pLU_dat(NULL),
	fIsFactorized(false),
	fNumThreads(0)	
{
	ExceptionT::GeneralFail("SPOOLESMatrixT_MT::SPOOLESMatrixT_MT", "not implemented");
}

/* destructor */
SPOOLESMatrixT_MT::~SPOOLESMatrixT_MT(void)
{
#ifndef OLD_CODE
	if (pLU_dat && !LU_MT_driver_free(&pLU_dat))
		ExceptionT::GeneralFail("SPOOLESMatrixT_MT::~SPOOLESMatrixT_MT",
			"error freeing SPOOLES data");
	pLU_dat = NULL;
#endif
}

/* clear values for next assembly */
void SPOOLESMatrixT_MT::Clear(void)
{
	const char caller[] = "SPOOLESMatrixT_MT::Clear";

	/* inherited */
	MSRMatrixT::Clear();

	/* check - must be serial */
	if (fTotNumEQ != fLocNumEQ)
		ExceptionT::GeneralFail(caller, "total equations (%d) != local equations (%d)",
			fTotNumEQ, fLocNumEQ);

#ifndef OLD_CODE
	/* delete existing data */
	if (pLU_dat && !LU_MT_driver_free(&pLU_dat))
		ExceptionT::GeneralFail(caller, "error freeing SPOOLES data");

	/* initialize new data */
	int matrix_type = SPOOLES_REAL;
	int symmetry_flag = (fSymmetric) ? SPOOLES_SYMMETRIC : SPOOLES_NONSYMMETRIC;
	int pivoting_flag = (fPivoting) ? SPOOLES_PIVOTING : SPOOLES_NO_PIVOTING;
	int seed = 1; /* any seed will do here */
	if (!LU_MT_driver_init(matrix_type, symmetry_flag, pivoting_flag, seed, fTotNumEQ, fNumThreads, &pLU_dat))
		ExceptionT::GeneralFail(caller, "error initializing SPOOLES data");
#endif /* OLD_CODE */

	/* reset flag */
	fIsFactorized = false;
}

/* assignment operator */
SPOOLESMatrixT_MT& SPOOLESMatrixT_MT::operator=(const SPOOLESMatrixT_MT&)
{
	const char caller[] = "SPOOLESMatrixT_MT::operator=";
	ExceptionT::GeneralFail(caller, "not implemented");
	return *this;
}

/** return a clone of self */
GlobalMatrixT* SPOOLESMatrixT_MT::Clone(void) const {
	return new SPOOLESMatrixT_MT(*this);
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* precondition matrix */
void SPOOLESMatrixT_MT::Factorize(void)
{
	const char caller[] = "SPOOLESMatrixT_MT::Factorize";

	/* quick exit */
	if (fIsFactorized) return;

	/* checks - must be serial */
	if (fTotNumEQ != fLocNumEQ)
		ExceptionT::GeneralFail(caller, "total equations (%d) != local equations (%d)",
			fTotNumEQ, fLocNumEQ);
	if (fStartEQ != 1)
		ExceptionT::GeneralFail(caller, "expecting first equation number to be 1 not %d",
			fStartEQ);

#ifndef OLD_CODE
	/* convert matrix to RCV */
	iArrayT r, c;
	dArrayT v;
	GenerateRCV(r, c, v, 1.0e-15);

	/* logging level:
		 0: nothing
		 1: scalar output (timing data) only
		>1: verbose */
	int msglvl = (fMessageLevel < 0) ? 0 : fMessageLevel; 

	/* compute factorization */
	int OK = LU_MT_driver_factorize(msglvl, SPOOLES_FILE,
		r.Length(), r.Pointer(), c.Pointer(), v.Pointer(), 
		&pLU_dat);

	if (OK != 1) ExceptionT::BadJacobianDet(caller, "LU driver returned %d", OK);

	/* set flag */
	fIsFactorized = true;
#endif
}

/* determine new search direction and put the results in result */
void SPOOLESMatrixT_MT::BackSubstitute(dArrayT& result)
{
	const char caller[] = "SPOOLESMatrixT_MT::BackSubstitute";
 
 	/* check */
	if (!fIsFactorized) ExceptionT::GeneralFail(caller, "matrix is not factorized");

	int msglvl = (fMessageLevel < 0) ? 0 : fMessageLevel; 
	 //  0: nothing
	 //  1: scalar output (timing data) only
	 // >1: verbose
	int OK = 1;

#ifdef OLD_CODE
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
		r.Length(), r.Pointer(), c.Pointer(), v.Pointer(), fNumThreads);

#else /* new driver */

	/* back substitute */
	OK = LU_MT_driver_solve(msglvl, SPOOLES_FILE, result.Pointer(), &pLU_dat);

#endif /* OLD_CODE */
 
	if (OK != 1) ExceptionT::BadJacobianDet(caller, "LU driver returned %d", OK);
}

#endif /* __SPOOLES_MT__ */
