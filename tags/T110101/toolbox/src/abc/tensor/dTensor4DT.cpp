/*
 * File: dTensor4DT.h 
 *
 */

/*
 * created      : PAK (05/25/97)
 * last modified: PAK (11/11/97)
 */

#include "dTensor4DT.h"

/*
 * Constructor
 */
dTensor4DT::dTensor4DT(void) { }
dTensor4DT::dTensor4DT(int dim0, int dim1, int dim2, int dim3):
	Tensor4DT<double>(dim0,dim1,dim2,dim3) { }
dTensor4DT::dTensor4DT(const dTensor4DT& source):
	Tensor4DT<double>(source) { }
