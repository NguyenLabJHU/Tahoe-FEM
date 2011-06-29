/*
 * File: dTensor3DT.cpp
 *
 */

/*
 * created      : PAK (05/24/97)
 * last modified: PAK (11/11/97)
 */

#include "dTensor3DT.h"

/*
 * Constructor
 */

using namespace Tahoe;

dTensor3DT::dTensor3DT(void) { }
dTensor3DT::dTensor3DT(int dim0, int dim1, int dim2):
	Tensor3DT<double>(dim0,dim1,dim2) { }
dTensor3DT::dTensor3DT(const dTensor3DT& source):
	Tensor3DT<double>(source) { }			

