/* $Id: SecantMethodT.cpp,v 1.6 2004-09-29 23:19:05 paklein Exp $ */
/* created: paklein (12/01/1998)*/
#include "SecantMethodT.h"
#include <math.h>

using namespace Tahoe;

/* constructor */
SecantMethodT::SecantMethodT(int max_iterations, double tolerance):
	fTol(tolerance),
	fMaxIts(max_iterations),
	ferr0(0.0),
	fcount(-1)
{
	/* check input values */
	if (fMaxIts < 1 || fTol < 0.0) ExceptionT::BadInputValue("SecantMethodT::SecantMethodT");
}

/* initialize the search with 2 intial guesses */
void SecantMethodT::Reset(double x1, double err1, double x2, double err2)
{
	/* initialize */
	Reset();
	NextPoint(x1, err1);
	NextPoint(x2, err2);
}

/* initialize */
void SecantMethodT::Reset(void)
{
	/* reset count */
	fcount = -1;
}

/* compute the next x guess, return 1 if converged or 0 if not */
double SecantMethodT::NextGuess(void) const
{
	const char caller[] = "SecantMethodT::NextGuess";
	if (fcount < 1) ExceptionT::GeneralFail(caller, "requires at least 2 initial guesses");

	double m = (ferr1 - ferr2)/(fx1 - fx2);
	if (fabs(m) < kSmall) ExceptionT::GeneralFail(caller, "bad slope = %g", m);

	double b = ferr2 - m*fx2;
	double x = -b/m;
	
	/* out of range -> contract */
//	if (x > fMaxStepSize || x < 0.0)
//		x = 0.5*(fx1 + fx2);
	
	return x;
}

/* try next point, returns 0 when converged */
SecantMethodT::StatusT SecantMethodT::NextPoint(double x, double err)
{
	if (fcount == -1) /* first guess */
	{
		fx1 = fx_best = x;
		ferr1 = ferr_best = err;
		fcount++;
		return kInit;
	}	
	else if (fcount == 0) /* second guess */
	{
		fx2 = x;
		ferr2 = err;
		ferr0 = (fabs(ferr1) > fabs(ferr2)) ? ferr2 : ferr1;
		if (fabs(ferr2) < fabs(ferr_best)) {
			fx_best = fx2;
			ferr_best = ferr2;
		}
		fcount++;
		return kInit;	
	}
	else
	{
		/* update best value */
		if (fabs(err) < fabs(ferr_best)) {
			ferr_best = err;
			fx_best = x;
		}

		/* stop */
		if (fabs(err) < fTol || fabs(err/ferr0) < fTol)
			return kConverged;
		else /* continue */
		{
			/* initialize flag */
			int give_up = 0;
	
			/* too many iterations */
			if (++fcount >= fMaxIts || fabs(fx1 - fx2) < kSmall)
				give_up = 1;
			else
			{
				if (fabs(ferr1) > fabs(err) && fabs(ferr1) > fabs(ferr2))
				{
					ferr1 = err;
					fx1 = x;
				}
				else if (fabs(ferr2) > fabs(err) && fabs(ferr2) > fabs(ferr1))
				{
					ferr2 = err;
					fx2 = x;
				}
				else
				{
					/* ferr1 and ferr2 don't bracket zero */
					if (ferr2*ferr1 > 0.0)
					{
						if (ferr1*err < 0)
						{
							ferr1 = err;
							fx1 = x;
						}
						else if (ferr2*err < 0)
						{
							ferr2 = err;
							fx2 = x;
						}
						else /* exit */
							give_up = 1;
					}
					else /* exit */
						give_up = 1;
				}
			}
			
			/* status */
			if (give_up == 1)
				return kFail;
			else
				return kContinue;
		}
	}
}
