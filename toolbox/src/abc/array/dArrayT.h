/* $Id: dArrayT.h,v 1.7 2003-11-21 22:41:30 paklein Exp $ */
/* created: paklein (08/11/1996) */
#ifndef _DARRAY_T_H_
#define _DARRAY_T_H_

/* base class */
#include "nArrayT.h"

namespace Tahoe {

class dArrayT: public nArrayT<double>
{
public:

	/** \name constructors */
	/*@{*/
	dArrayT(void);
	explicit dArrayT(int length);
	dArrayT(const dArrayT& source);

	/** construct an alias */
	dArrayT(int length, const double* p);
	/*@}*/

	/** \name assigment operators */
	/*@{*/
	dArrayT& operator=(const dArrayT& RHS);
	dArrayT& operator=(const double value);
	/*@}*/

	/** L2 norm of the vector */
	double Magnitude(void) const;

	/** \name create a unit vectors */
	/*@{*/
	dArrayT& UnitVector(const dArrayT& vector);

	/** scale this */
	dArrayT& UnitVector(void);
	/*@}*/
};

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

} // namespace Tahoe
 
#endif /* _DARRAY_T_H_ */
