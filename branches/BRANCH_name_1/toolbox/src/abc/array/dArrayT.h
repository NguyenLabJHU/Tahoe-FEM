/* $Id: dArrayT.h,v 1.4.4.2 2002-06-27 22:46:23 cjkimme Exp $ */
/* created: paklein (08/11/1996)                                          */

#ifndef _DARRAY_T_H_
#define _DARRAY_T_H_

/* base class */
#include "nArrayT.h"


namespace Tahoe {

class dArrayT: public nArrayT<double>
{
public:

	/* constructors */
	dArrayT(void);
	explicit dArrayT(int length);
	dArrayT(int length, double* p);
	dArrayT(const dArrayT& source);

	/* assigment operators */
	dArrayT& operator=(const dArrayT& RHS);
	dArrayT& operator=(const double value);

	/* L2 norm of the vector */
	double Magnitude(void) const;

	/* set this equal the unit vector in the direction of vector.  If no
	 * argument is passed in, *this is scaled */
	dArrayT& UnitVector(const dArrayT& vector);
	dArrayT& UnitVector(void);
};

} // namespace Tahoe

/* inlines */

/* assigment operators */
inline dArrayT& dArrayT::operator=(const dArrayT& RHS)
{
	nArrayT<double>::operator=(RHS);
	return *this;
}

inline dArrayT& dArrayT::operator=(const double value)
{
	nArrayT<double>::operator=(value);
	return *this;
}

inline dArrayT& dArrayT::UnitVector(const dArrayT& vector)
{
	SetToScaled(1.0/vector.Magnitude(), vector);
	return *this;
}

inline dArrayT& dArrayT::UnitVector(void)
{
	return UnitVector(*this);
}
 
#endif /* _DARRAY_T_H_ */
