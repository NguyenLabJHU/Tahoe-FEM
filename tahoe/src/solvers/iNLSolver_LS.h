/* $Id: iNLSolver_LS.h,v 1.2.2.3 2002-04-30 08:22:06 paklein Exp $ */
/* created: paklein (01/01/2001) */

#ifndef _I_NL_SOLVER_LS_H_
#define _I_NL_SOLVER_LS_H_

/* base class */
#include "NLSolver_LS.h"

/* direct members */
#include "dArray2DT.h"

/** nonlinear Newton solver with interactive console */
class iNLSolver_LS: public NLSolver_LS
{
public:

	/* constructor */
	iNLSolver_LS(FEManagerT& fe_manager, int group);

	/** solve the system over the current time increment */
	virtual int Solve(void);	
	
	/* execute commands */
	virtual bool iDoCommand(const CommandSpecT& command, StringT& line);

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
