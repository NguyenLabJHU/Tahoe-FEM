/* $Id: dTensor4DT.h,v 1.5 2002-07-13 00:07:00 cfoster Exp $ */
/* created paklein (05/25/97) */

#ifndef _D_TENSOR4D_T_H_
#define _D_TENSOR4D_T_H_

/* base class */
#include "Tensor4DT.h"
#include "dMatrixT.h"

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

	/* Converts 6x6 Matrix form of the Tangent Modulus to the 3x3x3x3
	 * form of the Tangent modulus */
	void ConvertTangentFrom2DTo4D(dTensor4DT& C, dMatrixT fc_ijkl);
	/* Converts 3x3x3x3 Tensor form of the Tangent Modulus 
	 * to the 6x6 matrix form of the Tangent modulus */
	/* Has not been tested !!!*/
	void ConvertTangentFrom4DTo2D(dTensor4DT& C, dMatrixT fc_ijkl);
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
