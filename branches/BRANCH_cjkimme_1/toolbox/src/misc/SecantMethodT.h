/* $Id: SecantMethodT.h,v 1.3 2002-09-12 16:40:20 paklein Exp $ */
/* created: paklein (12/01/1998)                                          */
/* SecantMethodT.h                                                        */

#ifndef _SECANT_METHOD_T_H_
#define _SECANT_METHOD_T_H_

#include "toolboxConstants.h"


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
};

} // namespace Tahoe 
#endif /* _SECANT_METHOD_T_H_ */
