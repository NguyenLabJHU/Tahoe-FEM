/**********************************************************************

	Error functions required to complete f2c libF77 support library
	
	IMT 4 Oct 97
	
	Modified to match changes in libF77; IMT 30 Nov 94
	
	No longer "necessary", but CW/MacOS erf() & erfc() behave badly 
	in the range (-0.5, 0.0).  More specifically, they return values 
	so inaccurate that I consider them out-right wrong in that range.
	
	All calls to the libf2c erf-related functions get mapped to
	my IMT_erf and IMT_erfc functions.
	
**********************************************************************/



#include <math.h>
#include "f2c.h"



/* Prototypes for local functions */

static void gcf( double *gammcf, double a, double x, double *gln );
static void gser( double *gamser, double a, double x, double *gln );
static double gammq( double a, double x );
static double gammp( double a, double x );
static double gammln( double xx );



double IMT_erf( double x )
{
	return  (x < 0.0) ? -gammp( 0.5, x*x ) : gammp( 0.5, x*x );
}

double IMT_erfc( double x )
{
	return  (x < 0.0) ? 1.0 + gammp( 0.5, x*x ) : gammq( 0.5, x*x );
}


/*
	Incomplete Gamma Function  P(a,x)
*/

static double gammp( double a, double x )
{
	double	gamser, gammcf, gln;

	if ( x < 0.0 || a <= 0.0 )
		/* This case can't happen when called by derf() or derfc() */ ;
		
	if (x < (a+1.0))
	{
		gser( &gamser, a, x, &gln );
		return  gamser;
	}
	else
		{
		gcf( &gammcf, a, x, &gln );
		return  1.0 - gammcf;
		}
}


/*
	Incomplete Gamma Function  Q(a,x)  =  1 - P(a,x)
*/

static double gammq( double a, double x )
{
	double 	gamser, gammcf, gln;

	if ( x < 0.0 || a <= 0.0 )
		/* This case can't happen when called by erf() or erfc() */ ;
		
	if ( x < (a + 1.0) )
	{
		gser( &gamser, a, x, &gln );
		return  1.0 - gamser;
	}
	else
	{
		gcf( &gammcf, a, x, &gln );
		return  gammcf;
	}
}


/*
	Returns following:
		gamser	-	Incomplete gamma function P(a,x) evaluated via series
		gln		-	ln(gamma function(a))
*/

#define  ITMAX 100
#define  EPS 3.0e-7

static void gser( double *gamser, double a, double x, double *gln )
{
	int		n;
	double	sum, del, ap;

	*gln = gammln(a);
	if (x <= 0.0)
	{
		if ( x < 0.0 )
			/* ERROR: This case can't happen when called by erf() or erfc() */ ;
		*gamser = 0.0;
		return;
	}
	else
	{
		ap = a;
		del = sum = 1.0/a;
		for ( n = 1; n <= ITMAX; n++ ) 
		{
			ap += 1.0;
			del *= x/ap;
			sum += del;
			if ( fabs(del) < fabs(sum)*EPS )
			{
				*gamser = sum*exp( -x + a*log(x) - *gln );
				return;
			}
		}
		/* ERROR: a too large, ITMAX too small.
		   This case can't happen when called by derf() or derfc() */ 
		return;
	}
}


/*
	Returns following:
		gamcf	-	Incomplete gamma function Q(a,x) evaluated via
					continued fraction
		gln		-	ln(gamma function(a))
*/

#define  ITMAX 100
#define  EPS 3.0e-7

static void gcf( double *gammcf, double a, double x, double *gln )
{
	int		n;
	double	gold = 0.0, g, fac = 1.0, b1 = 1.0,
			b0 = 0.0, anf, ana, an, a1, a0 = 1.0;

	*gln = gammln( a );
	a1 = x;
	for ( n = 1; n <= ITMAX; n++ )
	{
		an = (double) n;
		ana = an - a;
		a0 = (a1 + a0*ana)*fac;
		b0 = (b1 + b0*ana)*fac;
		anf = an*fac;
		a1 = x*a0 + anf*a1;
		b1 = x*b0 + anf*b1;
		if ( a1 )
		{
			fac = 1.0/a1;
			g = b1*fac;
			if ( fabs( (g-gold)/g ) < EPS )
			{
				*gammcf = exp( -x + a*log(x) - *gln )*g;
				return;
			}
			gold = g;
		}
	}
	/* ERROR: a too large, ITMAX too small.
	   This case can't happen when called by derf() or derfc() */ 
	return;
}

/*
	Natural log of gamma function.  Full accuracy obtained for xx > 1.
	For 0 < xx < 1, use reflection formula first.
*/

static double gammln( double xx )
{
	double		x, tmp, ser;
	int			j;

	static double	cof[6] = { 76.18009173, -86.505322033, 24.01409822,
							   -1.231739516, 0.120858003e-2, -0.536382e-5 };

	x = xx - 1.0;
	tmp = x + 5.5;
	tmp -= (x+0.5)*log(tmp);
	ser = 1.0;
	for (j = 0; j <= 5; j++)
		{
		x += 1.0;
		ser += cof[j]/x;
		}
	return  -tmp+log(2.50662827465*ser);
}




