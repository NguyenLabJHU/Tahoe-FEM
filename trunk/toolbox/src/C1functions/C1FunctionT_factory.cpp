/* $Id: C1FunctionT_factory.cpp,v 1.1 2004-01-27 19:07:23 paklein Exp $ */
#include "C1FunctionT.h"
#include <string.h>

/* subclasses supporting the factory method */
#include "PiecewiseLinearT.h"
#include "CubicSplineT.h"
#include "PowerLawT.h"

using namespace Tahoe;

/* factory method */
C1FunctionT* C1FunctionT::New(const char* name)
{
	if (strcmp(name, "piecewise_linear") == 0)
		return new PiecewiseLinearT;
	else if (strcmp(name, "cubic_spline") == 0)
		return new CubicSplineT;
	else if (strcmp(name, "power_law") == 0)
		return new PowerLawT;
	else
		return NULL;
}
