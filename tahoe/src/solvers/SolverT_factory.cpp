/* $Id: SolverT_factory.cpp,v 1.1.2.1 2004-02-05 18:47:17 paklein Exp $ */
#include "SolverT.h"
#include <string.h>

/* subclasses supporting the factory method */
#include "LinearSolver.h"
#include "NLSolver.h"
#include "PCGSolver_LS.h"

using namespace Tahoe;

/* factory method */
SolverT* SolverT::New(FEManagerT& fe_manager, const char* name)
{
	if (strcmp(name, "linear_solver") == 0)
		return new LinearSolver(fe_manager);
	else if (strcmp(name, "nonlinear_solver") == 0)
		return new NLSolver(fe_manager);
	else if (strcmp(name, "PCG_solver") == 0)
		return new PCGSolver_LS(fe_manager);
	else
		return NULL;
}
