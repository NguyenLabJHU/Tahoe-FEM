/* $Id: dTensor3DT.h,v 1.3 2002-07-05 22:26:21 paklein Exp $ */

/*
 * created      : PAK (05/24/97)
 * last modified: PAK (11/11/97)
 */

#ifndef _D_TENSOR3D_T_H_
#define _D_TENSOR3D_T_H_

/* base class */
#include "Tensor3DT.h"


namespace Tahoe {

class dTensor3DT: public Tensor3DT<double>
{
  public:

	/*
	 * Constructor
	 */
	dTensor3DT(void);
	dTensor3DT(int dim0, int dim1, int dim2);			
	dTensor3DT(const dTensor3DT& source);			

  	/*
  	 * Assignment operators
  	 */
  	dTensor3DT& operator=(const dTensor3DT& RHS);
  	dTensor3DT& operator=(const double value);

};

/* Inlines */

/*
 * Assignment operators
 */
inline dTensor3DT& dTensor3DT::operator=(const dTensor3DT& RHS)
{
	/* inherited */
	Tensor3DT<double>::operator=(RHS);
	return (*this);
}

inline dTensor3DT& dTensor3DT::operator=(const double value)
{
	/* inherited */
	Tensor3DT<double>::operator=(value);
	return (*this);
}

} // namespace Tahoe 
#endif /* _Z_TENSOR3D_T_H_ */
