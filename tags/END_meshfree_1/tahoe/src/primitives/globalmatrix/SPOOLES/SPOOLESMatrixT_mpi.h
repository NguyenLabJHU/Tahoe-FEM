/* $Id: SPOOLESMatrixT_mpi.h,v 1.2 2001-02-28 00:29:49 paklein Exp $ */
/* created: paklein (09/13/2000)                                          */

#ifndef _SPOOLES_MATRIX_T_MPI_H_
#define _SPOOLES_MATRIX_T_MPI_H_

/* base class */
#include "SPOOLESMatrixT.h"

/* library support options */
#ifdef __SPOOLES_MPI__
#ifdef __MPI__

class SPOOLESMatrixT_mpi: public SPOOLESMatrixT
{
public:

	/* constuctor */
	SPOOLESMatrixT_mpi(ostream& out, int check_code, bool symmetric,
		bool pivoting);

	/* assignment operator - not implemented */
	virtual GlobalMatrixT& operator=(const GlobalMatrixT& RHS);

protected:

	/* determine new search direction and put the results in result */
	virtual void BackSubstitute(dArrayT& result);
};
#endif /* __MPI__ */
#endif /*__SPOOLES_MPI__ */
#endif /* _SPOOLES_MATRIX_T_MPI_H_ */
