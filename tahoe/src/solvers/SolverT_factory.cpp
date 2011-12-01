/* $Id: SolverT_factory.cpp,v 1.3 2011-12-01 21:11:40 bcyansfn Exp $ */
#include "SolverT.h"
#include <cstring>

/* subclasses supporting the factory method */
#include "LinearSolver.h"
#include "NLSolver.h"
#include "PCGSolver_LS.h"
#include "NLSolver_LS.h"

using namespace Tahoe;

/* factory method */
SolverT* SolverT::New(FEManagerT& fe_manager, const char* name, int group)
{
	if (strcmp(name, "linear_solver") == 0)
		return new LinearSolver(fe_manager, group);
	else if (strcmp(name, "nonlinear_solver") == 0)
		return new NLSolver(fe_manager, group);
	else if (strcmp(name, "PCG_solver") == 0)
		return new PCGSolver_LS(fe_manager, group);
	else if (strcmp(name, "nonlinear_solver_LS") == 0)
		return new NLSolver_LS(fe_manager, group);
	else
		return NULL;
}
