/* $Id: SPOOLESMatrixT_mpi.h,v 1.3.4.1 2002-06-27 18:04:06 cjkimme Exp $ */
/* created: paklein (09/13/2000)                                          */

#ifndef _SPOOLES_MATRIX_T_MPI_H_
#define _SPOOLES_MATRIX_T_MPI_H_

/* base class */
#include "SPOOLESMatrixT.h"

/* library support options */
#ifdef __SPOOLES_MPI__
#ifdef __MPI__


namespace Tahoe {

class SPOOLESMatrixT_mpi: public SPOOLESMatrixT
{
public:

	/* constuctor */
	SPOOLESMatrixT_mpi(ostream& out, int check_code, bool symmetric,
		bool pivoting);

protected:

	/* determine new search direction and put the results in result */
	virtual void BackSubstitute(dArrayT& result);
};
} // namespace Tahoe 
#endif /* __MPI__ */
#endif /*__SPOOLES_MPI__ */
#endif /* _SPOOLES_MATRIX_T_MPI_H_ */
