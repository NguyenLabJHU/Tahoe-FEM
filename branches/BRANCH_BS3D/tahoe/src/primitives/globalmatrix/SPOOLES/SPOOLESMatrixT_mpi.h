/* $Id: SPOOLESMatrixT_mpi.h,v 1.6 2003-01-27 07:00:30 paklein Exp $ */
/* created: paklein (09/13/2000) */
#ifndef _SPOOLES_MATRIX_T_MPI_H_
#define _SPOOLES_MATRIX_T_MPI_H_

/* base class */
#include "SPOOLESMatrixT.h"

/* library support options */
#ifdef __SPOOLES_MPI__
#ifdef __TAHOE_MPI__

namespace Tahoe {

/* forward declarations */
class CommunicatorT;

class SPOOLESMatrixT_mpi: public SPOOLESMatrixT
{
public:

	/** constuctor */
	SPOOLESMatrixT_mpi(ostream& out, int check_code, bool symmetric,
		bool pivoting, CommunicatorT& comm);

protected:

	/* determine new search direction and put the results in result */
	virtual void BackSubstitute(dArrayT& result);
	
private:

	/** the MP layer */
	CommunicatorT& fComm;
};
} // namespace Tahoe 
#endif /* __TAHOE_MPI__ */
#endif /*__SPOOLES_MPI__ */
#endif /* _SPOOLES_MATRIX_T_MPI_H_ */
