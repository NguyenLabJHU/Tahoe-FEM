/* $Id: ExpCD_DRSolver.h,v 1.1.1.1 2001-01-29 08:20:34 paklein Exp $ */
/* created: paklein (08/19/1998)                                          */
/* quasi-static solver using dynamic relaxation with                      */
/* explicit central difference pseudo-dynamics                            */

#ifndef _EXPCD_DRSOLVER_H_
#define _EXPCD_DRSOLVER_H_

/* base class */
#include "SolverT.h"

/* direct members */
#include "fstreamT.h"

class ExpCD_DRSolver: public SolverT
{
public:

	/* constructor */
	ExpCD_DRSolver(FEManagerT& fe_manager);

	/* (re-)configure the global equation system */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);
	
	/* generate the solution for the current time sequence */
	 virtual void Run(void);

protected:

	/* advance to next load step. Returns 0 if there are no more
	 * steps. Overload to add class dependent initializations */
	virtual int Step(void);

	/* returns the appropriate iteration status flag for
	 * the given error measurement, based on the current
	 * iteration number, convergence tolerance, etc */
	int ExitIteration(double error);
	int ExitRelaxation(double error);

	/* form and solve the equation system - returns the magnitude of the
	 * residual */
	double SolveAndForm(void);
		// need this function ??????
	
	/* relax system - reform tangent at newtancount intervals */
	void Relax(int newtancount = 1);

private:

	/* set pseudo mass and damping parameters */
	void SetMass(void);
	double SetDamping(void);
	
protected:

	/* error management parameters */	
	int		fMaxIterations;
	double	fTolerance;		/* exponent - pow(10, tolerance);     */

	/* DR fudge parameters */
	double fMass_scaling;
	double fDamp_scaling;
	
	/* convergence history */
	int fOutputDOF; // <none> if {nodenum,dofnum} specifies
					// an inactive dof, or if nodenum < 1
	ofstreamT fhist_out;

	/* runtime error management data */
	double	fError0;
	
private:

	/* pseudo-dynamics vectors */
	dArrayT fDis;
	dArrayT fVel;
	dArrayT fAcc;
	
	dArrayT* pMass;
	
	/* time step */
	double fdt;
	
	/* "local" stiffness - based on initial config */
	dArrayT fK_0;
};

#endif /* _EXPCD_DRSOLVER_H_ */
