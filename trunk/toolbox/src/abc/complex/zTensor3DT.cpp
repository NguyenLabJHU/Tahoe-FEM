/* $Id: zTensor3DT.cpp,v 1.5 2002-10-20 22:38:52 paklein Exp $ */
/* created: PAK (05/19/97) */

#include "zTensor3DT.h"
#include "dTensor3DT.h"

/* constructor */

using namespace Tahoe;

zTensor3DT::zTensor3DT(void) { }
zTensor3DT::zTensor3DT(int dim0, int dim1, int dim2):
	Tensor3DT<ComplexT>(dim0,dim1,dim2) { }
zTensor3DT::zTensor3DT(const dTensor3DT& re, const dTensor3DT& im)
{
	toZ(re,im);	
}
zTensor3DT::zTensor3DT(const zTensor3DT& source): Tensor3DT<ComplexT>(source) { }

/*
 * Returning the Real and Imaginary parts
 */
void zTensor3DT::toRe(dTensor3DT& re) const
{
	/* dimension check */
	if (!SameDimensions(*this,re)) throw ExceptionT::kOutOfRange;

	/* ComplexT function */
	ComplexT::z_to_Re(*this, re);
}

void zTensor3DT::toIm(dTensor3DT& im) const
{
	/* dimension check */
	if (!SameDimensions(*this,im)) throw ExceptionT::kOutOfRange;

	/* ComplexT function */
	ComplexT::z_to_Im(*this, im);
}

zTensor3DT& zTensor3DT::toZ(const dTensor3DT& re, const dTensor3DT& im)
{
	/* dimension check */
	if (!SameDimensions(re,im)) throw ExceptionT::kOutOfRange;
	
	/* dimension */
	Dimension(re.Dim(0), re.Dim(1), re.Dim(2));

	/* ComplexT function */
	ComplexT::ReIm_to_z(re,im,*this);
	
	return (*this);
}
