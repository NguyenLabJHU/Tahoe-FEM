/* $Id: ModBessel.cpp,v 1.10 2003-11-21 22:41:27 paklein Exp $ */
/* created: dzeigle (4/18/2002) */

#include "ModBessel.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

/* constants */

using namespace Tahoe;

const int NUSE1 = 7;
const int NUSE2 = 8;
const double EPS = 1.0e-16;
const double FPMIN = 1.0e-30;
const int MAXIT = 10000;
const double XMIN = 2.0;
const double PI = 3.141592653589793;

/*
* Constructor
*/
ModBessel::ModBessel(double p)
{
	power = p;
}

/*
* Destructor
*/
ModBessel::~ModBessel()
{
}

/*
* Evaluation function
*/
double ModBessel::Eval(int d, double x) const
{
    double kn, dkn;
    
    kn = dkn = 0.0;
    
    bessk(x,power,&kn,&dkn);
    
    if (d==0)
        return kn;
    if (d==1)
        return dkn;
    else
    {
        cout << "\n**Error - noninteger derivative requested in ModBess::Eval.**\n";
        return 0.0;
    }
}
	 
/*
* I/O
*/
void ModBessel::Print(ostream& out) const
{
	/* parameters */
	out << " Scaling constant. . . . . . . . . . . . . . . . = " << power << '\n';
}

void ModBessel::PrintName(ostream& out) const
{
	out << "    Modified Bessel Function\n";
}

/*
* Returning values
*/
double ModBessel::Function(double x) const
{
	return Eval(0,x);
}

double ModBessel::DFunction(double x) const
{
	return Eval(1,x);
}

double ModBessel::DDFunction(double x) const
{
	cout << "\n Second derivative of the Bessel Function of the 3rd Kind not tabulated!\n";
	return 0.0*x;
}

/*
* Returning values in groups - derived classes should define
* their own non-virtual function called within this functon
* which maps in to out w/o requiring a virtual function call
* everytime. Default behavior is just to map the virtual functions
* above.
*/
dArrayT& ModBessel::MapFunction(const dArrayT& in,  dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;
	
	const double* pl = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = *pl++;	
		*pU++ = Eval(0,r);
	}
	return out;
}

dArrayT& ModBessel::MapDFunction(const dArrayT& in,  dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl  = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = *pl++;					
		*pdU++ = Eval(1,r);
	}
	return out;
}

dArrayT& ModBessel::MapDDFunction(const dArrayT& in,  dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl   = in.Pointer();
	double* pddU = out.Pointer();
	
	cout << "\n Second derivative of the Bessel Function of the 3rd Kind not tabulated!\n";
	for (int i = 0; i < in.Length(); i++)
	{
//		double r = *pl++;				
		*pl++;
		*pddU++ = 0.0;
	}
	return out;
}


// Evaluation of the Chebyshev function (see "Numerical Recipes in C", pg. 193.)
// Required for evaluation of Bessel functions of 1st and 2nd kind.
double ModBessel::chebev(double a, double b, double c[], int m, double x) const
{
    double d = 0.0, dd = 0.0, sv, y, y2;
    
    if ((x-a)*(x-b) > 0.0)
        cout << "\n**Error - x not in range in routine chebev.\n**";
    
    y2 = 2.0*(y=(2.0*x-a-b)/(b-a));
    for (int j=m-1; j>=1; j--)
    {
        sv = d;
        d = y2*d-dd+c[j];
        dd = sv;
    }
    
    return y*d-dd+0.5*c[0];
}

// Bessel functions of the 1st and 2nd kind. Taken from "Numerical Recipes in C", pg. 245.
// Analytic argument relies upon series formulation. Evalutes gamma function by
// Chebyshev expansion for |x| <= 0.5.
void ModBessel::beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi) const
{
    double xx;
    static double c1[] = {-1.142022680371168e0, 6.5165112670737e-3, 3.087090173086e-4, -3.4706269649e-6, 6.9437664e-9, 3.67795e-11, -1.356e-13};
    static double c2[] = {1.843740587300905e0, -7.68528408447867e-2, 1.2719271366546e-3, -4.9717367042e-6, -3.31261198e-8, 2.423096e-10, -1.702e-13, -1.49e-15};
    
    xx = 8.0*x*x-1.0;
    *gam1 = chebev(-1.0,1.0,c1,NUSE1,xx);
    *gam2 = chebev(-1.0,1.0,c2,NUSE2,xx);
    *gampl = *gam2-x*(*gam1);
    *gammi = *gam2+x*(*gam1);
}


// Returns the function value and derivative of the modified Bessel function.
// Borrowed portions from "Numerical Recipes in C", pg. 248. 
// The relative accuracy is within one or two significant digits of EPS. FPMIN is 
// a numerb close to the machine's smallest floating-point number.
void ModBessel::bessk(double x, double xnu, double *rk, double *rkp) const
{
    int i, l, nl;
    double a, a1, b, c, d, del, del1, delh, dels, e, fact, fact2, ff, gam1, gam2, gammi, gampl, h, p, pimu, q, q1, q2, qnew, ril, ripl, ritemp, rk1, rkmu, rktemp, s, sum, sum1, x2, xi, xi2, xmu, xmu2;
    
    if (x <= 0.0)
        cout << "\n**Error - bad arguments in bessk.\n**";
    if (xnu < 0.0)
        xnu = -xnu;
    
    nl = int(xnu+0.5);
    xmu = xnu-nl;
    xmu2 = xmu*xmu;
    xi = 1.0/x;
    xi2 = 2.0*xi;
    h = xnu*xi;
    if (h < FPMIN)
        h = FPMIN;
    b = xi2*xnu;
    d = 0.0;
    c = h;
    for (i=1; i<MAXIT; i++)
    {
        b += xi2;
        d = 1.0/(b+d);
        c = b+1.0/c;
        del = c*d;
        h = del*h;
        if (fabs(del-1.0)<EPS)
            break;
    }
    
    if (i>MAXIT)
        cout << "\n**Error - x too large in bessk; try asypmtotic expansion.\n**";
    ril = FPMIN;
    ripl = h*ril;
    fact = xnu*xi;
    for (l=nl; l>=1; l--)
    {
        ritemp = fact*ril+ripl;
        fact -= xi;
        ripl = fact*ritemp+ril;
        ril = ritemp;
    }
 
    if (x<XMIN)
    {
        x2 = 0.5*x;
        pimu  = PI*xmu;
        fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
        d = -log(x2);
        e = xmu*d;
        fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
        beschb(xmu, &gam1, &gam2, &gampl, &gammi);
        ff = fact*(gam1*cosh(e)+gam2*fact2*d);
        sum = ff;
        e = exp(e);
        p = 0.5*e/gampl;
        q = 0.5/(e*gammi);
        c = 1.0;
        d = x2*x2;
        sum1 = p;
        for (i=1; i<=MAXIT; i++)
        {
            ff = (i*ff+p+q)/(i*i-xmu2);
            c *= (d/i);
            p /= (i-xmu);
            q /= (i+xmu);
            del = c*ff;
            sum += del;
            del1 = c*(p-i*ff);
            sum1 += del1;
            if (fabs(del) < fabs(sum)*EPS)
                break;
        }
        if (i>MAXIT)
            cout << "\n**Error - bessk series failed to converge.\n**";
            rkmu = sum;
            rk1 = sum1*xi2;
    }
    else
    {
        b = 2.0*(1.0+x);
        d = 1.0/b;
        h = delh = d;
        q1 = 0.0;
        q2 = 1.0;
        a1 = 0.25-xmu2;
        q = c = a1;
        a = -a1;
        s = 1.0+q*delh;
        for (i=2; i<=MAXIT; i++)
        {
            a -= 2*(i-1);
            c = -a*c/i;
            qnew = (q1-b*q2)/a;
            q1 = q2;
            q2 = qnew;
            q += c*qnew;
            b += 2.0;
            d = 1.0/(b+a*d);
            delh = (b*d-1.0)*delh;
            h += delh;
            dels = q*delh;
            s += dels;
            if (fabs(dels/s) < EPS)
                break;
        }
        if (i>MAXIT)
            cout << "\n**Error - bessk: failure to converge in cf2.\n**";
        h = a1*h;
        rkmu = sqrt(PI/(2.0*x))*exp(-x)/s;
        rk1 = rkmu*(xmu+x+0.5-h)*xi;
    }

    for (i=1; i<=nl; i++)
    {
        rktemp = (xmu+i)*xi2*rk1+rkmu;
        rkmu = rk1;
        rk1 = rktemp;
    }
    *rk = rkmu;
    *rkp = xnu*xi*rkmu-rk1;
}


