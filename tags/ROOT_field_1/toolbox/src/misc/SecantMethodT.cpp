/* $Id: SecantMethodT.cpp,v 1.1.1.1 2001-01-25 20:56:25 paklein Exp $ */
/* created: paklein (12/01/1998)                                          */
/* SecantMethodT.cpp                                                      */

#include "SecantMethodT.h"

#include <math.h>

#include "ExceptionCodes.h"

/* constructor */
SecantMethodT::SecantMethodT(int max_iterations, double tolerance):
	fMaxIts(max_iterations),
	fTol(tolerance),
	ferr0(0.0),
	fcount(-1)
{
	/* check input values */
	if (fMaxIts < 1 || fTol < 0.0) throw eBadInputValue;	
}

/* initialize the search with 2 intial guesses */
void SecantMethodT::Reset(double x1, double err1, double x2, double err2)
{
	/* set values */
	fx1   = x1;
	ferr1 = err1;

	fx2   = x2;
	ferr2 = err2;
	
	/* set data */
	ferr0  = (fabs(ferr1) > fabs (ferr2)) ? fabs(ferr1) : fabs(ferr2);
	fcount = 0;	
}

/* compute the next x guess, return 1 if converged or 0 if not */
double SecantMethodT::NextGuess(void) const
{
	double m = (ferr1 - ferr2)/(fx1 - fx2);
	double b = ferr2 - m*fx2;

	return(-b/m);
}

/* try next point, returns 0 when converged */
int SecantMethodT::NextPoint(double x, double err)
{
	/* too many iterations */
	if (++fcount > fMaxIts) return -1;

	/* converged */
	if (fabs(err) < fTol || fabs(err/ferr0) < fTol)
		return 1;
	/* not converged */
	else
	{
		int give_up = 0;
	
		if ( fabs(ferr1) > fabs(err) && fabs(ferr1) > fabs(ferr2) )
		{
			ferr1   = err;
			fx1     = x;
		}
		else if ( fabs(ferr2) > fabs(err) && fabs(ferr2) > fabs(ferr1) )
		{
			ferr2   = err;
			fx2     = x;
		}
		else
		{
			/* ferr1 and ferr2 don't bracket zero */
			if (ferr2*ferr1 > 0.0)
			{
				if (ferr1*err < 0)
				{
					ferr1 = err;
					fx1   = x;
				}
				else if (ferr2*err < 0)
				{
					ferr2 = err;
					fx2   = x;
				}
				else /* exit */
					give_up = -1;
			}
			else /* exit */
				give_up = -1;
		}
		
		return give_up;
	}
}
