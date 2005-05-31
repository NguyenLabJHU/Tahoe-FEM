/* $Id: NLSolver_LSX.h,v 1.2.48.2 2005-05-31 06:08:55 paklein Exp $ */
#ifndef _NL_SOLVER_LSX_H_
#define _NL_SOLVER_LSX_H_

/* base class */
#include "NLSolver_LS.h"

namespace Tahoe {

/** temporary extension to NLSolver_LS which attempts to "push through"
 * trouble spots during the calculation */
class NLSolver_LSX: public NLSolver_LS
{
public:

	/** constructor */
	NLSolver_LSX(FEManagerT& fe_manager, int group);

	/** start solution step */
	virtual void InitStep(void);

protected:

	/** allows continuation of NLSolver_LSX::fMinStepRelError is less
	 * that NLSolver_LSX::fPuntTol */
	SolutionStatusT ExitIteration(double error, int iteration);

private:

	/** if min error less that this tolerance, move onto the
	 * next time step */
	double fPuntTol;

	/** double smallest error over the current solution step */
	double fMinStepRelError;
};

} /* namespace Tahoe */

#endif /* _NL_SOLVER_LSX_H_ */
