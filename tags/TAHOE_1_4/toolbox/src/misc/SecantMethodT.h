/* $Id: SecantMethodT.h,v 1.4 2004-05-26 09:32:10 paklein Exp $ */
/* created: paklein (12/01/1998) */
#ifndef _SECANT_METHOD_T_H_
#define _SECANT_METHOD_T_H_

/* direct members */
#include "dArray2DT.h"

namespace Tahoe {

class SecantMethodT
{
public:

	/* constructor */
	SecantMethodT(int max_iterations, double tolerance = 100*kSmall);

	/* initialize the search with 2 intial guesses */
	void Reset(double x1, double err1, double x2, double err2);

	/* return the next x guess */
	double NextGuess(void) const;
	
	/* try next point, returns 1 when converged, 0 if not converged,
	 * -1 on fail */
	int NextPoint(double x, double err);
	
	/** return the current number of iterations */
	int Iterations(void) const { return fcount; };

private:

	/* convergence tolerance */
	double fTol;
	int    fMaxIts;
	
	/* points */
	double fx1, ferr1;
	double fx2, ferr2;
	
	/* solution data */
	double ferr0;   // reference error
	int    fcount;  // number of iterations

	/** line search history */
	dArray2DT fSearchData;

	/** max step size */
	double fMaxStepSize;

	/** flag to run one more and then give up */
	bool fOneMore;
};

} /* namespace Tahoe */

#endif /* _SECANT_METHOD_T_H_ */
