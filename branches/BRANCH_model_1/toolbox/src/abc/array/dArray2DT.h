/* $Id: dArray2DT.h,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (07/16/1996)                                          */

#ifndef _DARRAY2D_T_H_
#define _DARRAY2D_T_H_

/* base class */
#include "nArray2DT.h"

/* forward declarations */
class dArrayT;

class dArray2DT: public nArray2DT<double>
{
public:

	/* constructors */
	dArray2DT(void);
	dArray2DT(int majordim, int minordim);
	dArray2DT(int majordim, int minordim, double* p);
	dArray2DT(const dArray2DT& source);

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

#endif /* _DARRAY2D_T_H */
