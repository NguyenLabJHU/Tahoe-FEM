/* $Id: nAutoFill2DT.h,v 1.2 2003-02-22 23:53:15 paklein Exp $ */
#ifndef _N_AUTO_FILL2D_T_H_
#define _N_AUTO_FILL2D_T_H_

/* base class */
#include "AutoFill2DT.h"

namespace Tahoe {

/** templated base class for arrays of number types. */
template <class nTYPE>
class nAutoFill2DT: public AutoFill2DT<nTYPE>
{
public:

	/** \name constructors */
	/*@{*/
	nAutoFill2DT(void) {};
	nAutoFill2DT(int majordim, int numchunks, int headroom, int maxminordim);
	/*@}*/

	/** \name operations with scalars */
	/*@{*/	
	nAutoFill2DT<nTYPE>& operator=(const nTYPE& value);  /**< set all elements in the array to value */
	nAutoFill2DT<nTYPE>& operator+=(const nTYPE& value); /**< add value to all elements */
	nAutoFill2DT<nTYPE>& operator-=(const nTYPE& value); /**< subtract value to all elements */
	nAutoFill2DT<nTYPE>& operator*=(const nTYPE& value); /**< multiply all elements by value */
	nAutoFill2DT<nTYPE>& operator/=(const nTYPE& value); /**< divide all elements by value */
	/*@}*/	
};

/*************************************************************************
 * Implementation
 *************************************************************************/

template <class nTYPE>
nAutoFill2DT<nTYPE>::nAutoFill2DT(int majordim, int numchunks, int headroom, int maxminordim):
	AutoFill2DT<nTYPE>(majordim, numchunks, headroom, maxminordim)
{

}

template <class nTYPE>
nAutoFill2DT<nTYPE>& nAutoFill2DT<nTYPE>::operator=(const nTYPE& value)
{
	int md = MajorDim();
	for (int i = 0; i < md; i++)
	{
		int dim = MinorDim(i);
		nTYPE* p = (*this)(i);
		for (int j = 0; j < dim; j++)
			*p++ = value;
	}
	return *this;
}

template <class nTYPE>
nAutoFill2DT<nTYPE>& nAutoFill2DT<nTYPE>::operator+=(const nTYPE& value)
{
	int md = MajorDim();
	for (int i = 0; i < md; i++)
	{
		int dim = MinorDim(i);
		nTYPE* p = (*this)(i);
		for (int j = 0; j < dim; j++)
			*p++ += value;
	}
	return *this;
}

template <class nTYPE>
nAutoFill2DT<nTYPE>& nAutoFill2DT<nTYPE>::operator-=(const nTYPE& value)
{
	int md = MajorDim();
	for (int i = 0; i < md; i++)
	{
		int dim = MinorDim(i);
		nTYPE* p = (*this)(i);
		for (int j = 0; j < dim; j++)
			*p++ -= value;
	}
	return *this;
}

template <class nTYPE>
nAutoFill2DT<nTYPE>& nAutoFill2DT<nTYPE>::operator*=(const nTYPE& value)
{
	int md = MajorDim();
	for (int i = 0; i < md; i++)
	{
		int dim = MinorDim(i);
		nTYPE* p = (*this)(i);
		for (int j = 0; j < dim; j++)
			*p++ *= value;
	}
	return *this;
}

template <class nTYPE>
nAutoFill2DT<nTYPE>& nAutoFill2DT<nTYPE>::operator/=(const nTYPE& value)
{
	int md = MajorDim();
	for (int i = 0; i < md; i++)
	{
		int dim = MinorDim(i);
		nTYPE* p = (*this)(i);
		for (int j = 0; j < dim; j++)
			*p++ /= value;
	}
	return *this;
}

} /* namespace Tahoe */

#endif /* _N_AUTO_FILL2D_T_H_ */
