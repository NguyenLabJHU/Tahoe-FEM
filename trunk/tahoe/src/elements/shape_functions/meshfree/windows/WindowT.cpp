/* $Id: WindowT.cpp,v 1.3 2004-10-12 00:20:26 paklein Exp $ */
#include "WindowT.h"
#include "dArrayT.h"
#include "dArray2DT.h"

using namespace Tahoe;

/* coverage test */
int WindowT::Covers(const dArray2DT& x_n, const dArrayT& x, const dArray2DT& param_n, ArrayT<bool>& covers) const
{
	dArrayT x_n_i;
	dArrayT param_n_i;
	int count = 0;
	for (int i = 0; i < x_n.MajorDim(); i++)
	{
		x_n.RowAlias(i, x_n_i);
		param_n.RowAlias(i, param_n_i);
		covers[i] = Covers(x_n_i, x, param_n_i);
		if (covers[i])
			count++;
	}

	return count;
}

/* compute spherical support size as batch */
void WindowT::SphericalSupportSize(const dArray2DT& param_n, ArrayT<double>& support_size) const
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
void WindowT::RectangularSupportSize(const dArray2DT& param_n, dArray2DT& support_size) const
{
#if __option(extended_errorcheck)
	if (support_size.MajorDim() != param_n.MajorDim())
		ExceptionT::SizeMismatch("WindowT::SphericalSupportSize");
#endif

	dArrayT param, support;
	for (int i = 0; i < support_size.MajorDim(); i++) {
		param_n.RowAlias(i, param);
		support_size.RowAlias(i, support);
		RectangularSupportSize(param, support);
	}	
}
