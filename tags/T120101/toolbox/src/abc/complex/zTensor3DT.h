/*
 * File: zTensor3DT.h 
 *
 */

/*
 * created      : PAK (05/19/97)
 * last modified: PAK (11/11/97)
 */

#ifndef _Z_TENSOR3D_T_H_
#define _Z_TENSOR3D_T_H_

/* base class */
#include "Tensor3DT.h"

/* direct members */
#include "ComplexT.h"

/* forward declarations */
class dTensor3DT;

class zTensor3DT: public Tensor3DT<ComplexT>
{
  public:

	/*
	 * Constructor
	 */
	zTensor3DT(void);
	zTensor3DT(int dim0, int dim1, int dim2);
	zTensor3DT(const dTensor3DT& re, const dTensor3DT& im);
	zTensor3DT(const zTensor3DT& source);

 	/*
	 * Assigment operators
	 */
	zTensor3DT& operator=(const zTensor3DT& RHS);
	zTensor3DT& operator=(const ComplexT& value);
			
	/*
  	 * Returning the Real and Imaginary parts
  	 */
  	void toRe(dTensor3DT& re) const;
  	void toIm(dTensor3DT& im) const;
  	zTensor3DT& toZ(const dTensor3DT& re, const dTensor3DT& im);	
};

/*
 * Assigment operators
 */
inline zTensor3DT& zTensor3DT::operator=(const zTensor3DT& RHS)
{
	Tensor3DT<ComplexT>::operator=(RHS);
	return(*this);
}

inline zTensor3DT& zTensor3DT::operator=(const ComplexT& value)
{
	Tensor3DT<ComplexT>::operator=(value);
	return(*this);
}

#endif /* _Z_TENSOR3D_T_H_ */
