/* $Id: dTensor4DT.h,v 1.4 2002-07-05 22:26:21 paklein Exp $ */
/* created paklein (05/25/97) */

#ifndef _D_TENSOR4D_T_H_
#define _D_TENSOR4D_T_H_

/* base class */
#include "Tensor4DT.h"

namespace Tahoe {

/** fourth order tensor class for double's. Most functionality
 * is inherited from the base class Tensor4DT. */
class dTensor4DT: public Tensor4DT<double>
{
  public:

	/** \name constructors */
	/*@{*/
	dTensor4DT(void);
	dTensor4DT(int dim0, int dim1, int dim2, int dim3);			
	dTensor4DT(const dTensor4DT& source);			
	/*@}*/

  	/** \name assignment operators */
	/*@{*/  	
	dTensor4DT& operator=(const dTensor4DT& RHS);
  	dTensor4DT& operator=(const double value);
	/*@}*/
};

inline dTensor4DT& dTensor4DT::operator=(const dTensor4DT& RHS)
{
	/* inherited */
	Tensor4DT<double>::operator=(RHS);
	return *this;
}

inline dTensor4DT& dTensor4DT::operator=(const double value)
{
	/* inherited */
	Tensor4DT<double>::operator=(value);
	return *this;
}

} // namespace Tahoe 
#endif /* _D_TENSOR4D_T_H_ */
