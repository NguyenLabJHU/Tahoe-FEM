/* $Id: GaussianWindowT.cpp,v 1.1 2001-06-13 20:53:56 paklein Exp $ */

#include "GaussianWindowT.h"
#include "ExceptionCodes.h"

/* constructor */
GaussianWindowT::GaussianWindowT(double dilation_scaling, double sharpening_factor):
	fDilationScaling(dilation_scaling),
	fSharpeningFudgeFactor(sharpening_factor)
{
	if (fDilationScaling < 0.0 || fSharpeningFudgeFactor < 0.0)
		throw eBadInputValue;
}
