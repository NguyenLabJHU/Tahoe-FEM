/* $Id: SecantMethodT.cpp,v 1.5 2004-05-26 09:32:10 paklein Exp $ */
/* created: paklein (12/01/1998)*/
#include "SecantMethodT.h"
#include <math.h>

using namespace Tahoe;

/* constructor */
SecantMethodT::SecantMethodT(int max_iterations, double tolerance):
	fTol(tolerance),
	fMaxIts(max_iterations),
	ferr0(0.0),
	fcount(-1),
	fSearchData(fMaxIts, 2),
	fMaxStepSize(2.5),
	fOneMore(false)
{
	/* check input values */
	if (fMaxIts < 1 || fTol < 0.0) throw ExceptionT::kBadInputValue;	
}

/* initialize the search with 2 intial guesses */
void SecantMethodT::Reset(double x1, double err1, double x2, double err2)
{
	/* clear history */
	fSearchData = 0.0;
	fOneMore = false;	

	/* set values */
	fx1 = x1;
	ferr1 = err1;

	fx2 = x2;
	ferr2 = err2;
	
	/* set data */
	ferr0 = (fabs(ferr1) > fabs(ferr2)) ? ferr2 : ferr1;
	fcount = 1;	

	/* store history */
	fSearchData(0,0) = fx1;
	fSearchData(0,1) = ferr1;

	fSearchData(1,0) = fx2;
	fSearchData(1,1) = ferr2;
}

/* compute the next x guess, return 1 if converged or 0 if not */
double SecantMethodT::NextGuess(void) const
{
	double m = (ferr1 - ferr2)/(fx1 - fx2);
	double b = ferr2 - m*fx2;
	double x = -b/m;
	
	/* out of range -> contract */
	if (x > fMaxStepSize || x < 0.0)
		x = 0.5*(fx1 + fx2);
	
	return x;
}

/* try next point, returns 0 when converged */
int SecantMethodT::NextPoint(double x, double err)
{
	/* stop */
	if (fOneMore) return 1;

	/* stop */
	if (fabs(err) < fTol || fabs(err/ferr0) < fTol)
		return 1;
	else /* continue */
	{
		/* initialize flag */
		int give_up = 0;
	
		/* too many iterations */
		if (++fcount >= fMaxIts || fabs(fx1 - fx2) < kSmall)
			give_up = 1;
		else
		{
			/* store */
			fSearchData(fcount, 0) = x;
			fSearchData(fcount, 1) = err;

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
	
		/* look for best so far */
		if (give_up)
		{
			/* find "best" step */
			int best = 0;
			if (fabs(fSearchData(best,0)) < kSmall) best = 1; // skip zero step
			double err_best = fabs(fSearchData(best,1));
			for (int i = best + 1; i < fSearchData.MajorDim(); i++)
			{
				double err_test = fabs(fSearchData(i,1));
				if (err_test < err_best)
				{
					err_best = err_test;
					best = i;
				}
			}

			/* set flag */
			fOneMore = true;
		}
		
		/* at least one more round */
		return 0;
	}
}
