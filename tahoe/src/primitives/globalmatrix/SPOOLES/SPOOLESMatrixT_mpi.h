/* $Id: SPOOLESMatrixT_mpi.h,v 1.10 2005-04-13 21:50:13 paklein Exp $ */
/* created: paklein (09/13/2000) */
#ifndef _SPOOLES_MATRIX_T_MPI_H_
#define _SPOOLES_MATRIX_T_MPI_H_

/* base class */
#include "SPOOLESMatrixT.h"

/* library support options */
#ifdef __SPOOLES_MPI__
#ifdef __TAHOE_MPI__

namespace Tahoe {

class SPOOLESMatrixT_mpi: public SPOOLESMatrixT
{
public:

	/** constuctor */
	SPOOLESMatrixT_mpi(ostream& out, int check_code, bool symmetric,
		bool pivoting, int message_level, const CommunicatorT& comm);

	/** destructor */
	~SPOOLESMatrixT_mpi(void);

	/** clear values for next assembly */
	virtual void Clear(void);

	/** assignment operator */
	SPOOLESMatrixT_mpi& operator=(const SPOOLESMatrixT_mpi& rhs);

protected:

	/** precondition matrix */
	virtual void Factorize(void);

	/** determine new search direction and put the results in result */
	virtual void BackSubstitute(dArrayT& result);

};
} // namespace Tahoe 
#endif /* __TAHOE_MPI__ */
#endif /*__SPOOLES_MPI__ */
#endif /* _SPOOLES_MATRIX_T_MPI_H_ */
