/* $Id: NLSolver_LS.h,v 1.7 2004-01-05 07:07:19 paklein Exp $ */
/* created: paklein (08/18/1999) */

#ifndef _NL_SOLVER_LS_H_
#define _NL_SOLVER_LS_H_

/* base class */
#include "NLSolver.h"

/* direct members */
#include "dArray2DT.h"

namespace Tahoe {

/** nonlinear Newton solver with line search */
class NLSolver_LS: public NLSolver
{
public:

	/* constructor */
	NLSolver_LS(FEManagerT& fe_manager, int group);

	/* form and solve the equation system - returns the magnitude of the
	 * residual */
	virtual double SolveAndForm(int& iteration);

	/* console */
	virtual bool iDoVariable(const StringT& variable, StringT& line);

protected:

	/* apply system update (socket for line searching) */
	virtual void Update(const dArrayT& update, const dArrayT* residual);

private:

	/* return the line search weight function for the given step size.
	 * The degrees of freedom:
	 *
	 *	G = R(d_i+1).delta_d */
	 double GValue(double step);

private:

	/* line search parameters */
	int    fSearchIterations;
	double fOrthogTolerance;
	double fMaxStepSize;

	/* work space */
	dArrayT fR; // store first residual

	/* line search data */
	dArrayT   fUpdate;     // full update vector
	double    s_current;   // current step size
	dArray2DT fSearchData; // line search history
};

} // namespace Tahoe 
#endif /* _NL_SOLVER_LS_H_ */
