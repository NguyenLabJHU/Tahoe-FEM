/* $Id: NLK0Solver.h,v 1.1.1.1 2001-01-29 08:20:34 paklein Exp $ */
/* created: paklein (10/01/1996)                                          */
/* Solver that forms and factorizes the tangent matrix only once          */
/* (for now), with line search updates                                    */

#ifndef _NL_K0_SOLVER_H_
#define _NL_K0_SOLVER_H_

/* base class */
#include "NLSolver.h"

/* direct members */
#include "CCSMatrixT.h"

class NLK0Solver: public NLSolver
{
public:

	/* constructor */
	NLK0Solver(FEManagerT& fe_manager);

	/* generate the solution for the current time sequence */
	 virtual void RunK0(void); //not used

protected:

	/* form and solve - returns the magnitude of the residual */
	virtual double SolveAndForm(bool junk);

private:

	/* perform linesearching to find the best step length */
	void InitLineSearch(void);
	void TermLineSearch(void);
	void LineSearchUpdate(double s_upper);
	
	/* return the line search weight function for the given step size.
	 * The degrees of freedom:
	 *
	 *	G = R(d_i+1).delta_d */
	 double GValue(double& step);

private:

	/* casted pointer to global matrix */
	CCSMatrixT*	pCCSLHS;
	
	/* last positive definite tangent (factorized) */
	int			fFormTangent; //no switching back, yet
	CCSMatrixT	fLastTangent;

	/* full update vector */
	dArrayT		fUpdate;

	/* line search data */
	double		s_current;			/* current step size */
	int			fRecursionDepth;	/* depth of GValue recursion */

};

#endif /* _NL_K0_SOLVER_H_ */
