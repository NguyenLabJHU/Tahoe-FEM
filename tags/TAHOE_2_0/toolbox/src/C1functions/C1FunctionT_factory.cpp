/* $Id: C1FunctionT_factory.cpp,v 1.8 2004-06-29 16:15:43 paklein Exp $ */
#include "C1FunctionT.h"
#include <string.h>

/* subclasses supporting the factory method */
#include "PiecewiseLinearT.h"
#include "CubicSplineT.h"
#include "PowerLawT.h"
#include "LinearT.h"
#include "LennardJones612.h"
#include "SmithFerrante.h"
#include "ModSmithFerrante.h"
#include "LinearExponentialT.h"
#include "GaoJi.h"
#include "GaoJi2.h"
#include "GaoVicky.h"
#include "SF2.h"
#include "CosineT.h"

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
	else if (strcmp(name, "modified_Smith-Ferrante") == 0)
		return new ModSmithFerrante;
	else if (strcmp(name, "linear_exponential") == 0)
		return new LinearExponentialT;
	else if (strcmp(name, "Gao-Ji") == 0)
		return new GaoJi;
	else if (strcmp(name, "Gao-Ji_2") == 0)
		return new GaoJi2;
	else if (strcmp(name, "Gao-Nguyen") == 0)
		return new GaoVicky;
	else if (strcmp(name, "Smith-Ferrante_2") == 0)
		return new SF2;
	else if (strcmp(name, "cosine") == 0)
		return new CosineT;
	else
		return NULL;
}
