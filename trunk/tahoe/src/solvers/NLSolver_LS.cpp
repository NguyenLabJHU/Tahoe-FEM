/* $Id: NLSolver_LS.cpp,v 1.15 2004-09-09 23:54:55 paklein Exp $ */
/* created: paklein (08/18/1999) */
#include "NLSolver_LS.h"

#include <iostream.h>
#include <math.h>

#include "toolboxConstants.h"
#include "ExceptionT.h"
#include "FEManagerT.h"

using namespace Tahoe;

/* constructor */
NLSolver_LS::NLSolver_LS(FEManagerT& fe_manager, int group):
	NLSolver(fe_manager, group),
	fSearchIterations(3),
	fOrthogTolerance(0.25),
	fMaxStepSize(2.5)
{
	SetName("nonlinear_solver_LS");

	/* set console */
	iAddVariable("line_search_iterations", fSearchIterations);
	iAddVariable("line_search_tolerance", fOrthogTolerance);
	iAddVariable("max_step_size", fMaxStepSize);
}

/* do one iteration of the solution procedure */
void NLSolver_LS::Iterate(void)
{	
	/* store residual */
	fR = fRHS;

	/* solve equation system */
	if (!fLHS->Solve(fRHS)) ExceptionT::BadJacobianDet("NLSolver_LS::Iterate");

	/* apply update to system - using line search */
	fRHS_lock = kOpen;
	Update(fRHS, &fR);									
}

/* console */
bool NLSolver_LS::iDoVariable(const StringT& variable, StringT& line)
{
	/* inherited */
	bool result = NLSolver::iDoVariable(variable, line);
	if (result)
	{
		/* need to reallocate */
		if (variable == "line_search_iterations")
			fSearchData.Dimension(fSearchIterations, 2);
	}
	return result;
}

/* describe the parameters needed by the interface */
void NLSolver_LS::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	NLSolver::DefineParameters(list);

	/* line search iterations */
	ParameterT line_search_iterations(fSearchIterations, "line_search_iterations");
	line_search_iterations.SetDefault(fSearchIterations);
	list.AddParameter(line_search_iterations);

	/* line search orthogonality tolerance */
	ParameterT line_search_tolerance(fOrthogTolerance, "line_search_tolerance");
	line_search_tolerance.SetDefault(fOrthogTolerance);
	list.AddParameter(line_search_tolerance);

	/* maximum step size */
	ParameterT max_step(fMaxStepSize, "max_step");
	max_step.SetDefault(fMaxStepSize);
	list.AddParameter(max_step);
}

/* accept parameter list */
void NLSolver_LS::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	NLSolver::TakeParameterList(list);

	/* extract line search parameters */
	fSearchIterations = list.GetParameter("line_search_iterations");
	fOrthogTolerance = list.GetParameter("line_search_tolerance");
	fMaxStepSize = list.GetParameter("max_step");

	/* allocate space for history */
	fSearchData.Dimension(fSearchIterations, 2);
	fSearchData = 0.0;
}

/*************************************************************************
* Protected
*************************************************************************/

void NLSolver_LS::Update(const dArrayT& update, const dArrayT* residual)
{
	/* inherited (no line search) */
	if (fSearchIterations == 0)
	{
		NLSolver::Update(update, residual);
		return;
	}

	/* initialize */
	fUpdate     = update;
	s_current   = 0.0;
	fSearchData = 0.0;

	/* first secant point */
	double s_a;
	double G_a;
	if (residual)
	{
		s_a = 0.0;
		G_a = InnerProduct(fUpdate, *residual);
	}
	else
	{
		s_a = 0.5;
		G_a = GValue(s_a);
	}	
	fSearchData(0,0) = s_a;
	fSearchData(0,1) = G_a;

	/* start with check full step */
	double s_b = 1.0;
	double G_b = GValue(s_b);
#if 1
	int cuts = 0;
	while (cuts++ < 10 && fabs(G_a) > kSmall && fabs(G_b/G_a) > 1.0e6) /* too big */
	{
		s_b = 0.5*(s_b + s_a);
		G_b = GValue(s_b);
	}
	cout << " init:" << setw(2) << cuts - 1;
	if (cuts == 10) throw ExceptionT::kBadJacobianDet;
#endif
	fSearchData(1,0) = s_b;
	fSearchData(1,1) = G_b;

	/* minimize G with secant method */		
	double G_0 = (fabs(G_a) > fabs(G_b)) ? G_b : G_a;
	int count = 2;
	double G_new;
	bool give_up = false;
	do {
		double m = (G_a - G_b)/(s_a - s_b);
		double b = G_b - m*s_b;
		double s_new = -b/m;
		
		/* out of range -> contract */
		if (s_new > fMaxStepSize || s_new < 0.0)
			s_new = 0.5*(s_a + s_b);

		/* update and compute test value */				
		G_new = GValue(s_new);
		fSearchData(count, 0) = s_new;
		fSearchData(count, 1) = G_new;

		/* secant search */
		if (fabs(G_a) > fabs(G_new) && fabs(G_a) > fabs(G_b))
		{
			G_a = G_new;
			s_a = s_new;
			give_up = false;
		}
		else if (fabs(G_b) > fabs(G_new) && fabs(G_b) > fabs(G_a))
		{
			G_b = G_new;
			s_b = s_new;
			give_up = false;
		}
		else
		{
			/* G_a and G_b don't bracket zero */
			if (G_b*G_a > 0)
			{
				if (G_a*G_new < 0)
				{
					G_a = G_new;
					s_a = s_new;
				}
				else if (G_b*G_new < 0)
				{
					G_b = G_new;
					s_b = s_new;
				}
				else /* exit */
					give_up = true;
			}
			else /* exit */
				give_up = true;
		
		}
		
		/* max iterations or s_a == s_b */
		if (++count >= fSearchIterations || fabs(s_a - s_b) < kSmall) 
			give_up = true;
		
	} while (fabs(G_new) > fZeroTolerance &&
	         fabs(G_new/G_0) > fOrthogTolerance &&
			!give_up);

	/* best step on fail */
	if (give_up)
	{
		/* find "best" step */
		int best = 0;
		if (fabs(fSearchData(best,0)) < kSmall) best = 1; // skip zero step
		double G_best = fabs(fSearchData(best,1));
		for (int i = best + 1; i < count; i++)
		{
			double G_test = fabs(fSearchData(i,1));
			if (G_test < G_best)
			{
				G_best = G_test;
				best = i;
			}
		}
	
		/* set to "best" */
		GValue(fSearchData(best,0));
	}
	
	/* no good */
	if (fabs(s_current) < kSmall) throw ExceptionT::kBadJacobianDet;

	/* write results */
	cout << " LS: " << count << setw(kDoubleWidth) << s_current << " | ";
} 	

/*************************************************************************
* Private
*************************************************************************/

/* return the line search weight function for the given step size */
double NLSolver_LS::GValue(double step)
{
	/* scale update vector */
	fRHS = fUpdate;

	/* scale accounting for current step update */
	fRHS *= (step - s_current);
	s_current = step;
	
	/* compute residual */
	fFEManager.Update(Group(), fRHS);
	fRHS = 0.0;
	try { fFEManager.FormRHS(Group()); }

	catch (ExceptionT::CodeT error)
	{
		cout << "\n NLSolver_LS::GValue: caught exception: " << error << endl;
		throw error;
	}
	
	return InnerProduct(fUpdate, fRHS);
}
