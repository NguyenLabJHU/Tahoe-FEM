/* $Id: SPOOLESMatrixT.cpp,v 1.21 2004-03-16 06:56:35 paklein Exp $ */
/* created: paklein (09/13/2000) */
#include "SPOOLESMatrixT.h"

/* library support options */
#ifdef __SPOOLES__

#include "ElementMatrixT.h"
#include "SPOOLES.h"

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
	fPivoting(pivoting)
{

}

SPOOLESMatrixT::SPOOLESMatrixT(const SPOOLESMatrixT& source):
	MSRMatrixT(source)
{
	ExceptionT::GeneralFail("SPOOLESMatrixT::SPOOLESMatrixT", "not implemented");
}

/* assignment operator */
GlobalMatrixT& SPOOLESMatrixT::operator=(const SPOOLESMatrixT& rhs)
{
#pragma unused(rhs)

	ExceptionT::GeneralFail("SPOOLESMatrixT::operator=", "not implemented");
	return *this;
}

/* assignment operator */
GlobalMatrixT& SPOOLESMatrixT::operator=(const GlobalMatrixT& rhs)
{
	const char caller[] = "SPOOLESMatrixT::operator=";

#ifdef __NO_RTTI__
	ExceptionT::GeneralFail(caller, "requires RTTI");
#endif

	const SPOOLESMatrixT* sp = TB_DYNAMIC_CAST(const SPOOLESMatrixT*, &rhs);
	if (!sp)  ExceptionT::GeneralFail(caller, "cast const SPOOLESMatrixT* failed");
	return operator=(*sp);
}

/** return a clone of self */
GlobalMatrixT* SPOOLESMatrixT::Clone(void) const
{
	SPOOLESMatrixT* new_mat = new SPOOLESMatrixT(*this);
	return new_mat;
}

/*************************************************************************
 * Protected
 *************************************************************************/
	
/* determine new search direction and put the results in result */
void SPOOLESMatrixT::BackSubstitute(dArrayT& result)
{
	const char caller[] = "SPOOLESMatrixT::BackSubstitute";

	/* check */
	if (fTotNumEQ != fLocNumEQ)
		ExceptionT::GeneralFail(caller, "total equations (%d) != local equations (%d)",
			fTotNumEQ, fLocNumEQ);

	/* convert matrix to RCV */
	iArrayT r, c;
	dArrayT v;
	GenerateRCV(r, c, v, 1.0e-15);

 
	 int msglvl = 0; //  0: nothing
	 //  1: scalar output (timing data) only
	 // >1: verbose
	int matrix_type = SPOOLES_REAL;
	int symmetry_flag = (fSymmetric) ? SPOOLES_SYMMETRIC : SPOOLES_NONSYMMETRIC;
	int pivoting_flag = (fPivoting) ? SPOOLES_PIVOTING : SPOOLES_NO_PIVOTING;
	int seed = 1; /* any seed will do here */
	int OK;
 
#ifdef __SPOOLES_MT__
	/* call MT driver, __NUM_THREADS__ comes from make macro */
	OK = LU_MT_driver(msglvl, SPOOLES_FILE, matrix_type, symmetry_flag,
		pivoting_flag, seed, result.Length(), result.Pointer(),
		r.Length(), r.Pointer(), c.Pointer(), v.Pointer(), __NUM_THREADS__);
#else /* not __SPOOLES_MT__ */
	OK = LU_serial_driver(msglvl, SPOOLES_FILE, matrix_type, symmetry_flag,
		pivoting_flag, seed, result.Length(), result.Pointer(),
		r.Length(), r.Pointer(), c.Pointer(), v.Pointer());
#endif /* __SPOOLES_MT__ */
 
	if (OK != 1) ExceptionT::BadJacobianDet(caller, "LU driver returned %d", OK);
}

#endif /* __SPOOLES__ */
