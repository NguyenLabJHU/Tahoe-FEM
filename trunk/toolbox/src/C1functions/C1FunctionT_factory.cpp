/* $Id: C1FunctionT_factory.cpp,v 1.4 2004-03-24 01:56:34 paklein Exp $ */
#include "C1FunctionT.h"
#include <string.h>

/* subclasses supporting the factory method */
#include "PiecewiseLinearT.h"
#include "CubicSplineT.h"
#include "PowerLawT.h"
#include "LinearT.h"
#include "LennardJones612.h"
#include "SmithFerrante.h"

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
	else if (strcmp(name, "linear_function") == 0)
		return new LinearT;
	else if (strcmp(name, "Lennard-Jones_6-12") == 0)
		return new LennardJones612;
	else if (strcmp(name, "Smith-Ferrante") == 0)
		return new SmithFerrante;
	else
		return NULL;
}
