/*
 * File: dTensor4DT.h 
 *
 */

/*
 * created      : PAK (05/25/97)
 * last modified: PAK (11/11/97)
 */

#ifndef _D_TENSOR4D_T_H_
#define _D_TENSOR4D_T_H_

/* base class */
#include "Tensor4DT.h"

class dTensor4DT: public Tensor4DT<double>
{
  public:

	/*
	 * Constructor
	 */
	dTensor4DT(void);
	dTensor4DT(int dim0, int dim1, int dim2, int dim3);			
	dTensor4DT(const dTensor4DT& source);			

  	/*
  	 * Assignment operators
  	 */
  	dTensor4DT& operator=(const dTensor4DT& RHS);
  	dTensor4DT& operator=(const double value);

};

/* Inlines */

/*
 * Assignment operators
 */
inline dTensor4DT& dTensor4DT::operator=(const dTensor4DT& RHS)
{
	/* inherited */
	Tensor4DT<double>::operator=(RHS);
	return (*this);
}

inline dTensor4DT& dTensor4DT::operator=(const double value)
{
	/* inherited */
	Tensor4DT<double>::operator=(value);
	return (*this);
}

#endif /* _Z_TENSOR4D_T_H_ */
