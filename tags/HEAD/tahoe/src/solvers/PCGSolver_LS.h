/* $Id: PCGSolver_LS.h,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (08/19/1999)                                          */

#ifndef _PCG_SOLVER_LS_H_
#define _PCG_SOLVER_LS_H_

/* base class */
#include "NLSolver.h"

/* direct members */
#include "dArray2DT.h"

class PCGSolver_LS: public NLSolver
{
public:

	/* constructor */
	PCGSolver_LS(FEManagerT& fe_manager);

	/* (re-)configure the global equation system */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);
	
protected:

	/* apply system update (socket for line searching) */
	virtual void Update(const dArrayT& update,
		const dArrayT* residual);

	/* form and solve the equation system - returns the magnitude
	 * of the residual */
	virtual double SolveAndForm(bool newtangent);

private:

	/* find new search conjugate search direction */
	void CGSearch(void);

	/* return the line search weight function for the given step size.
	 * The degrees of freedom:
	 *
	 *	G = R(d_i+1).delta_d */
	double GValue(double step);
	
private:

	/* conjugate gradient parameters */
	int fRestart;

	/* line search parameters */
	int    fSearchIterations;
	double fOrthogTolerance;
	double fMaxStepSize;

	/* work space */
	dArrayT fR;      // residual
	dArrayT fR_last; // last residual
	dArrayT fu_last; // last update
	dArrayT fdiff_R; // residual difference

	//TEMP
	int fPreconditioner; //flag to reform preconditioner

	/* line search data */
	dArrayT   fUpdate;     // full update vector
	double    s_current;   // current step size
	dArray2DT fSearchData; // line search history
};

#endif /* _PCG_SOLVER_LS_H_ */
