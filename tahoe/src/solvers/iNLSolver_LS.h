/* $Id: iNLSolver_LS.h,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (01/01/2001)                                          */

#ifndef _I_NL_SOLVER_LS_H_
#define _I_NL_SOLVER_LS_H_

/* base class */
#include "NLSolver_LS.h"

/* direct members */
#include "dArray2DT.h"

class iNLSolver_LS: public NLSolver_LS
{
public:

	/* constructor */
	iNLSolver_LS(FEManagerT& fe_manager);

	/* generate the solution for the current time sequence */
	virtual void Run(void);
	
	/* execute commands */
	virtual bool iDoCommand(const StringT& command, StringT& line);

private:

	/* apply flags */
	virtual void Update(const dArrayT& update, const dArrayT* residual);

	/* commands */
	bool DoStep(int max_steps);
	bool DoInitStep(void);
	IterationStatusT DoIterate(int max_count);
	
private:

	/* ON/OFF flags */
	bool fFormTangent;	
	bool fLineSearch;	

	/* solution states */
	IterationStatusT fIterationStatus;
};

#endif /* _I_NL_SOLVER_LS_H_ */
