/* $Id: DRSolver.h,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: PAK/CBH (10/03/1996)                                          */

#ifndef _DRSOLVER_H_
#define _DRSOLVER_H_

/* base class */
#include "NLSolver.h"

/* forward declarations */
class CCSMatrixT;

class DRSolver: public NLSolver
{
public:

	/* constructor */
	DRSolver(FEManagerT& fe_manager);

	/* configure the global equation system */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);
	
	/* generate the solution for the current time sequence */
	 virtual void Run(void);

private:

	/* compute the pseudo-mass */
	void ComputeMass(void);
	void ComputeVelocity(void);
	void ComputeDamping(void);

private:

	dArrayT fMass;
	dArrayT fVel;
	dArrayT fDisp;
	double  fDamp;
	
	/* work space */
	dArrayT fKd;		

	int fNumEquations;
	CCSMatrixT* fCCSLHS; //dynamically casted base class data
	 	  	
};

#endif /* _DRSOLVER_H_ */
