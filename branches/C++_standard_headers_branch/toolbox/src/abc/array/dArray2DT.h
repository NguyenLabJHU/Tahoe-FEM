/* $Id: dArray2DT.h,v 1.3 2003-11-21 22:41:30 paklein Exp $ */
/* created: paklein (07/16/1996) */

#ifndef _DARRAY2D_T_H_
#define _DARRAY2D_T_H_

/* base class */
#include "nArray2DT.h"

namespace Tahoe {

/* forward declarations */
class dArrayT;

class dArray2DT: public nArray2DT<double>
{
public:

	/** \name constructors */
	/*@{*/
	dArray2DT(void);
	dArray2DT(int majordim, int minordim);
	dArray2DT(const dArray2DT& source);

	/** construct an alias */
	dArray2DT(int majordim, int minordim, const double* p);
	/*@}*/

	/* assignment operators */
	dArray2DT& operator=(const dArray2DT& RHS);
	dArray2DT& operator=(const double value);
};

/* inlines */

/* assigment operators */
inline dArray2DT& dArray2DT::operator=(const dArray2DT& RHS)
{
	nArray2DT<double>::operator=(RHS);
	return *this;
}

inline dArray2DT& dArray2DT::operator=(const double value)
	{
	nArray2DT<double>::operator=(value);
	return *this;
}

} // namespace Tahoe 
#endif /* _DARRAY2D_T_H */
