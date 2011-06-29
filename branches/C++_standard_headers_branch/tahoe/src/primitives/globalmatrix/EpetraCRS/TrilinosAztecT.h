/* $Id: TrilinosAztecT.h,v 1.1 2006-11-01 05:15:32 paklein Exp $ */
#ifndef _TRILINOS_AZTEC_T_H_
#define _TRILINOS_AZTEC_T_H_

/* library support */
#ifdef __TRILINOS__

/* base class */
#include "EpetraCRSMatrixT.h"

namespace Tahoe {

/** interface Aztec solver in Trilinos */
class TrilinosAztecT: public EpetraCRSMatrixT
{
public:

	/** constructor */
	TrilinosAztecT(ostream& out, int check_code, const CommunicatorT& comm);

	/** copy constructor */
	TrilinosAztecT(const TrilinosAztecT& rhs);

	/** assignment operator */
	TrilinosAztecT& operator=(const TrilinosAztecT& rhs);

	/** return a clone of self */
	virtual GlobalMatrixT* Clone(void) const;

protected:

	/** solution driver  */
	virtual void BackSubstitute(dArrayT& result);

};

} /* namespace Tahoe */

#endif /* __TRILINOS__ */
#endif /* _TRILINOS_AZTEC_T_H_ */
