/* $Id: NLSolver_LS.h,v 1.9 2004-09-09 23:54:55 paklein Exp $ */
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

	/** constructor */
	NLSolver_LS(FEManagerT& fe_manager, int group);

	/** do one iteration of the solution procedure */
	virtual void Iterate(void);

	/* console */
	virtual bool iDoVariable(const StringT& variable, StringT& line);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

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

	/** \name line search parameters */
	/*@{*/
	int    fSearchIterations;
	double fOrthogTolerance;
	double fMaxStepSize;
	/*@}*/

	/* work space */
	dArrayT fR; // store first residual

	/** \name line search data */
	/*@{*/
	dArrayT   fUpdate;     /**< full update vector */
	double    s_current;   /**< current step size */
	dArray2DT fSearchData; /**< line search history */
	/*@}*/
};

} // namespace Tahoe 
#endif /* _NL_SOLVER_LS_H_ */
