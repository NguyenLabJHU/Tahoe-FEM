/* $Id: WindowGPT.cpp,v 1.1 2004-06-29 17:15:15 raregue Exp $ */
#include "WindowGPT.h"
#include "dArrayT.h"
#include "dArray2DT.h"

using namespace Tahoe;

/* compute spherical support size as batch */
void WindowGPT::SphericalSupportSize(const dArray2DT& param_n, ArrayT<double>& support_size) const
{
#if __option(extended_errorcheck)
	if (support_size.Length() != param_n.MajorDim())
		ExceptionT::SizeMismatch("WindowT::SphericalSupportSize");
#endif

	dArrayT param;
	for (int i = 0; i < support_size.Length(); i++) {
		param_n.RowAlias(i, param);
		support_size[i] = SphericalSupportSize(param);
	}	
}

/* compute rectangular support size as batch */
void WindowGPT::RectangularSupportSize(const dArray2DT& param_n, dArray2DT& support_size) const
{
#if __option(extended_errorcheck)
	if (support_size.MajorDim() != param_n.MajorDim())
		ExceptionT::SizeMismatch("WindowT::SphericalSupportSize");
#endif

	dArrayT param, support;
	for (int i = 0; i < support_size.Length(); i++) {
		param_n.RowAlias(i, param);
		support_size.SetRow(i, RectangularSupportSize(param));
	}	
}
