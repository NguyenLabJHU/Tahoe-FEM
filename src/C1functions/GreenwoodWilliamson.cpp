/* $Id: GreenwoodWilliamson.cpp,v 1.4 2002-02-04 17:38:28 dzeigle Exp $ */

#include "GreenwoodWilliamson.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "ModBessel.h"
#include "ConHyperGeom.h"
#include "GenHyperGeom.h"

/* constants */
const double PI = 2.0*acos(0.0);

/* functions */
static double dom1int(double s, double h, double mean, double stdev)
{
	double xval=pow((s-mean)/stdev,2.0);
	double phi=(exp(-0.5*xval))/(stdev*sqrt(2.0*PI));
	
	return (pow((s-h),1.5)*phi);
}

static double ddom1intdh(double s, double h, double mean, double stdev)
{
	double xval=pow((s-mean)/stdev,2.0);
	double phi=(exp(-0.5*xval))/(stdev*sqrt(2.0*PI));
	
	return (-1.5*sqrt(s-h)*phi);
}

/*
* constructors
*/
GreenwoodWilliamson::GreenwoodWilliamson(double MU, double SIGMA):
	fM(MU), fS(SIGMA) { }

/*
* I/O
*/
void GreenwoodWilliamson::Print(ostream& out) const
{
	/* parameters */
	out << " Potential parameters:\n";
	out << "      MU = " << fM << '\n';
	out << "      SIGMA = " << fS << '\n';
}

void GreenwoodWilliamson::PrintName(ostream& out) const
{
	out << "    Greenwood and Williamson\n";
}

/*
* Returning values
*/
double GreenwoodWilliamson::Function(double x) const
{
	cout << "*** ERROR! Greenwood and Williamson potential unavailable.";
	return (-1.0e10*x);
}

double GreenwoodWilliamson::DFunction(double x) const
{
	double xval, diff, yval;
	
	if (fS==0)
	{
		cout << "\n*** Bad SIGMA value in GreenwoodWilliamson.cpp.\n";
		throw eBadInputValue;
	}
	else
	{
		diff = x - fM;
		xval = pow(diff/(2.0*fS),2.0);
		
		if (diff > 0)
		{
			ModBessel k1(0.25), k3(0.75), k5(1.25);
			double f0 = sqrt(diff/PI)*exp(-xval)/(8.0*fS);
			double f1 = 2.0*(pow(diff,2.0) + 2.0*pow(fS,2.0))*k1.Function(xval);
			double f2 = -pow(diff,2.0)*(k3.Function(xval) + k5.Function(xval));
	
			yval = f0*(f1 + f2);
		}
		else if (diff == 0)
			yval = 0.906402*pow(fS,1.5)/(pow(2.0,0.25)*sqrt(PI)); 
		else
		{	/* decompose "negative domain": [h,mu] + [mu,inf]	*/
			/*													*/
			/* no closed form solution; use Gaussian quadrature */
			/* on the function:									*/
			/*			(s-h)^(1.5) phi(s)						*/
			double dom1=0.0;
			double a = x;
			double b = fM;
			double x1, x2, x3, x4;
			double A1, A2, A3, A4;
			double temp1, temp2, temp3, temp4, temp5, temp6;
			double temp7, temp8, temp9, temp10, temp11;
			double ap2, ap3, ap4;
			double bp2, bp3;
			double x1p2, x1p3;
			double x2p2, x2p3;
			double x3p2;
			double x4p2, x4p3;
			
			x1 = (35.0*a-sqrt(525.0+70.0*sqrt(30.0))*(a-b)+35.0*b)/70.0;
			x2 = (35.0*a+sqrt(525.0+70.0*sqrt(30.0))*(a-b)+35.0*b)/70.0;
			x3 = (35.0*a-sqrt(525.0-70.0*sqrt(30.0))*(a-b)+35.0*b)/70.0;
			x4 = (35.0*a+sqrt(525.0-70.0*sqrt(30.0))*(a-b)+35.0*b)/70.0;
			
			ap2 = pow(a,2.0);
			ap3 = pow(a,3.0);
			ap4 = pow(a,4.0);
			bp2 = pow(b,2.0);
			bp3 = pow(b,3.0);
			x1p2 = pow(x1,2.0);
			x1p3 = pow(x1,3.0);
			x2p2 = pow(x2,2.0);
			x2p3 = pow(x2,3.0);
			x3p2 = pow(x3,2.0);
			x4p2 = pow(x4,2.0);
			x4p3 = pow(x4,3.0);
			
			temp1 = 12.0*a*x2*(x2-x3)*(x2-x4)*(x3-x4)*x4;
			temp2 = 4.0*ap3*(2.0*x2p3-2.0*x2*x3p2+x4*x3p2-x4p3);
			temp3 = ap4*(-6.0*x2p2+6.0*x2*x3+3.0*x4*(x4-x3));
			temp4 = 6.0*ap2*x3*(-2.0*x2p3+2.0*x2p2*x3+x4p2*(x4-x3));
			temp5 = -12.0*x2*(x2-x3)*(x2-x4)*(x3-x4)*x4;
			temp6 = 3*bp3*(2.0*x2p2-2.0*x2*x3+(x3-x4)*x4);
			temp7 = 6.0*b*x3*(2.0*x2p3-2.0*x2p2*x3+(x3-x4)*x4p2);
			temp8 = bp2*(-8.0*x2p3+8.0*x2*x3p2-4.0*x3p2*x4+4.0*x4p3);
			temp9 = (x2-x3)*(x2-x4)*(x3-x4)*x4;
			temp10 = x1p2*(2.0*x2p2-2.0*x2*x3+(x3-x4)*x4);
			temp11 = x1*(-2.0*x2p2*x3-x3p2*x4+x4p3+x2*(2.0*x3p2+x3*x4-x4p2));
			A1 = (temp1+temp2+temp3+temp4+b*(temp5+temp6+temp7+temp8))/(12.0*(x1-x2)*(temp9+temp10+temp11));
			
			temp1 = -12.0*a*x1*(x1-x3)*(x1-x4)*(x3-x4)*x4;
			temp2 = 3.0*ap4*(2.0*x1p2-2.0*x1*x3+(x3-x4)*x4);
			temp3 = 6.0*ap2*x3*(2.0*x1p3-2.0*x1p2*x3+(x3-x4)*x4p2);
			temp4 = ap3*(-8.0*x1p3+8.0*x1*x3p2-4.0*x3p2*x4+4.0*x4p3);
			temp5 = 12.0*x1*(x1-x3)*(x1-x4)*(x3-x4)*x4;
			temp6 = 4*bp2*(2.0*x1p3-2.0*x1*x3p2+x3p2*x4-x4p3);
			temp7 = bp3*(-6.0*x1p2+6.0*x1*x3+3.0*x4*(x4-x3));
			temp8 = 6.0*b*x3*(-2.0*x1p3+2.0*x1p2*x3+x4p2*(x4-x3));
			temp9 = (x2-x3)*(x2-x4)*(x3-x4)*x4;
			temp10 = x1p2*(2.0*x2p2-2.0*x2*x3+(x3-x4)*x4);
			temp11 = x1*(-2.0*x2p2*x3-x3p2*x4+x4p3+x2*(2.0*x3p2+x3*x4-x4p2));
			A2 = (temp1+temp2+temp3+temp4+b*(temp5+temp6+temp7+temp8))/(12.0*(x1-x2)*(temp9+temp10+temp11));
			
			temp1 = 12.0*a*x1*x2*(x1-x4)*(x2-x4)*x4;
			temp2 = -3.0*ap4*(2.0*x1*x2-x1*x4-x2*x4+x4p2);
			temp3 = 4.0*ap3*(x1p2*(2.0*x2-x4)+x1*x2*(2.0*x2-x4)-x2p2*x4+x4p3);
			temp4 = -6*ap2*(x1*x4p2*(x4-x2)+x2*x4p2*(x4-x2)+x1p2*(2.0*x2p2-x4p2));
			temp5 = -12.0*x1*x2*(x1-x4)*(x2-x4)*x4;
			temp6 = 3.0*bp3*(2.0*x1*x2-x1*x4-x2*x4+x4p2);
			temp7 = -4.0*bp2*(x1p2*(2.0*x2-x4)+x1*x2*(2.0*x2-x4)-x2p2*x4+x4p3);
			temp8 = 6.0*b*(x1*x4p2*(x4-x2)+x2*x4p2*(x4-x2)+x1p2*(2.0*x2p2-x4p2));
			temp9 = (x2-x3)*(x2-x4)*(x3-x4)*x4;
			temp10 = x1p2*(2.0*x2p2-2.0*x2*x3+x4*(x3-x4));
			temp11 = x1*(-2.0*x2p2*x3-x3p2*x4+x4p3+x2*(2.0*x3p2+x3*x4-x4p2));
			A3 = (temp1+temp2+temp3+temp4+b*(temp5+temp6+temp7+temp8))/(12.0*x3*(temp9+temp10+temp11));
			
			temp1 = -12.0*a*x1*x2*(x1-x3)*(x2-x3)-3.0*ap4*(x1+x2-x3);
			temp2 = -6.0*x3*ap2*(x1p2+x1*(x2-x3)+x2*(x2-x3));
			temp3 = 4.0*ap3*(x1p2+x1*x2+x2p2-x3p2);
			temp4 = 12.0*x1*x2*(x1-x3)*(x2-x3)+3.0*bp3*(x1+x2-x3);
			temp5 = 6.0*b*x3*(x1p2+x1*x2+x2p2-x1*x3-x2*x3);
			temp6 = -4.0*bp2*(x1p2+x1*x2+x2p2-x3p2);
			temp7 = (x2-x3)*(x2-x4)*(x3-x4)*x4;
			temp8 = x1p2*(2.0*x2p2-2.0*x2*x3+(x3-x4)*x4);
			temp9 = x1*(-2.0*x2p2*x3-x3p2*x4+x4p3+x2*(2.0*x3p2+x3*x4-x4p2));
			A4 = (temp1+temp2+temp3+b*(temp4+temp5+temp6))/(12.0*(temp7+temp8+temp9));
			
			dom1 = A1*dom1int(x1,x,fM,fS)+A2*dom1int(x2,x,fM,fS)+A3*dom1int(x3,x,fM,fS)+A4*dom1int(x4,x,fM,fS);

			double dom2 = 0.0;
			ConHyperGeom fCHG3(-0.75,0.5), fCHG1(-0.25,1.5);
			GenHyperGeom fGHG(0.5,1.0,1.75,2.25);
			
			temp1 = sqrt(fM-x)/(5.0*fS*sqrt(2.0*PI));
			temp2 = 5.0*pow(2.0,0.25)*pow(fS,1.5)/sqrt(fM-x);
			temp3 = fS*0.906402*fCHG3.Function(-2.0*xval);
			temp4 = sqrt(2.0)*(fM-x)*0.919063*fCHG1.Function(-2.0*xval);
			temp5 = -2.0*pow((fM-x),2.0)*fGHG.Function(-2.0*xval);
			dom2 = temp1*(temp2*(temp3+temp4)+temp5);

			yval = dom1+dom2;
		}
	}
		return yval;
}

double GreenwoodWilliamson::DDFunction(double x) const
{
	double xval, diff, yval;
	
	if (fS==0)
	{
		cout << "\n*** Bad SIGMA value in GreenwoodWilliamson.cpp.\n";
		throw eBadInputValue;
	}
	else
	{
		diff = x - fM;
		xval = pow(diff/(2.0*fS),2.0);
		
		if (diff > 0)
		{
			ModBessel k1(0.25), k3(0.75), k5(1.25), kn1(-0.25), k7(1.75), kn3(-0.75), k9(2.25);
			double f0 = sqrt(fS*(x-fM))*exp(-xval);
			double f1 = 2.0*(pow(x-fM,2.0) + 2.0*pow(fS,2.0))*k1.Function(xval);
			double f2 = -pow(x-fM,2.0)*(k3.Function(xval) + k5.Function(xval));
			
			double df0 = -exp(-xval)*(pow(x-fM,2.0)-pow(fS,2.0))/(2.0*pow(fS,1.5)*sqrt(x-fM));
			double df1 = (x-fM)*(4.0*k1.Function(xval)-(2.0*xval+1.0)*(kn3.Function(xval)+k5.Function(xval)));
			double df2sum = kn1.Function(xval)+k1.Function(xval)+k7.Function(xval)+k9.Function(xval);
			double df2 = (x-fM)*(-2.0*(k3.Function(xval)+k5.Function(xval))+xval*df2sum);
	
			yval = df0*(f1+f2)+f0*(df1+df2);
		}
		else if (diff == 0)
			yval = -3.0*1.22542*sqrt(fS/PI)/pow(2.0,1.75);
		else
		{	/* decompose "negative domain": [h,mu] + [mu,inf]	*/
			/*													*/
			/* no closed form solution; use Gaussian quadrature */
			/* on the function:									*/
			/*			(s-h)^(1.5) phi(s)						*/
			double dom1=0.0;
			double a = x;
			double b = fM;
			double x1, x2, x3, x4;
			double A1, A2, A3, A4;
			double temp1, temp2, temp3, temp4, temp5, temp6;
			double temp7, temp8, temp9, temp10, temp11;
			double ap2, ap3, ap4;
			double bp2, bp3;
			double x1p2, x1p3;
			double x2p2, x2p3;
			double x3p2;
			double x4p2, x4p3;
			
			x1 = (35.0*a-sqrt(525.0+70.0*sqrt(30.0))*(a-b)+35.0*b)/70.0;
			x2 = (35.0*a+sqrt(525.0+70.0*sqrt(30.0))*(a-b)+35.0*b)/70.0;
			x3 = (35.0*a-sqrt(525.0-70.0*sqrt(30.0))*(a-b)+35.0*b)/70.0;
			x4 = (35.0*a+sqrt(525.0-70.0*sqrt(30.0))*(a-b)+35.0*b)/70.0;
			
			ap2 = pow(a,2.0);
			ap3 = pow(a,3.0);
			ap4 = pow(a,4.0);
			bp2 = pow(b,2.0);
			bp3 = pow(b,3.0);
			x1p2 = pow(x1,2.0);
			x1p3 = pow(x1,3.0);
			x2p2 = pow(x2,2.0);
			x2p3 = pow(x2,3.0);
			x3p2 = pow(x3,2.0);
			x4p2 = pow(x4,2.0);
			x4p3 = pow(x4,3.0);
			
			temp1 = 12.0*a*x2*(x2-x3)*(x2-x4)*(x3-x4)*x4;
			temp2 = 4.0*ap3*(2.0*x2p3-2.0*x2*x3p2+x4*x3p2-x4p3);
			temp3 = ap4*(-6.0*x2p2+6.0*x2*x3+3.0*x4*(x4-x3));
			temp4 = 6.0*ap2*x3*(-2.0*x2p3+2.0*x2p2*x3+x4p2*(x4-x3));
			temp5 = -12.0*x2*(x2-x3)*(x2-x4)*(x3-x4)*x4;
			temp6 = 3*bp3*(2.0*x2p2-2.0*x2*x3+(x3-x4)*x4);
			temp7 = 6.0*b*x3*(2.0*x2p3-2.0*x2p2*x3+(x3-x4)*x4p2);
			temp8 = bp2*(-8.0*x2p3+8.0*x2*x3p2-4.0*x3p2*x4+4.0*x4p3);
			temp9 = (x2-x3)*(x2-x4)*(x3-x4)*x4;
			temp10 = x1p2*(2.0*x2p2-2.0*x2*x3+(x3-x4)*x4);
			temp11 = x1*(-2.0*x2p2*x3-x3p2*x4+x4p3+x2*(2.0*x3p2+x3*x4-x4p2));
			A1 = (temp1+temp2+temp3+temp4+b*(temp5+temp6+temp7+temp8))/(12.0*(x1-x2)*(temp9+temp10+temp11));
			
			temp1 = -12.0*a*x1*(x1-x3)*(x1-x4)*(x3-x4)*x4;
			temp2 = 3.0*ap4*(2.0*x1p2-2.0*x1*x3+(x3-x4)*x4);
			temp3 = 6.0*ap2*x3*(2.0*x1p3-2.0*x1p2*x3+(x3-x4)*x4p2);
			temp4 = ap3*(-8.0*x1p3+8.0*x1*x3p2-4.0*x3p2*x4+4.0*x4p3);
			temp5 = 12.0*x1*(x1-x3)*(x1-x4)*(x3-x4)*x4;
			temp6 = 4*bp2*(2.0*x1p3-2.0*x1*x3p2+x3p2*x4-x4p3);
			temp7 = bp3*(-6.0*x1p2+6.0*x1*x3+3.0*x4*(x4-x3));
			temp8 = 6.0*b*x3*(-2.0*x1p3+2.0*x1p2*x3+x4p2*(x4-x3));
			temp9 = (x2-x3)*(x2-x4)*(x3-x4)*x4;
			temp10 = x1p2*(2.0*x2p2-2.0*x2*x3+(x3-x4)*x4);
			temp11 = x1*(-2.0*x2p2*x3-x3p2*x4+x4p3+x2*(2.0*x3p2+x3*x4-x4p2));
			A2 = (temp1+temp2+temp3+temp4+b*(temp5+temp6+temp7+temp8))/(12.0*(x1-x2)*(temp9+temp10+temp11));
			
			temp1 = 12.0*a*x1*x2*(x1-x4)*(x2-x4)*x4;
			temp2 = -3.0*ap4*(2.0*x1*x2-x1*x4-x2*x4+x4p2);
			temp3 = 4.0*ap3*(x1p2*(2.0*x2-x4)+x1*x2*(2.0*x2-x4)-x2p2*x4+x4p3);
			temp4 = -6*ap2*(x1*x4p2*(x4-x2)+x2*x4p2*(x4-x2)+x1p2*(2.0*x2p2-x4p2));
			temp5 = -12.0*x1*x2*(x1-x4)*(x2-x4)*x4;
			temp6 = 3.0*bp3*(2.0*x1*x2-x1*x4-x2*x4+x4p2);
			temp7 = -4.0*bp2*(x1p2*(2.0*x2-x4)+x1*x2*(2.0*x2-x4)-x2p2*x4+x4p3);
			temp8 = 6.0*b*(x1*x4p2*(x4-x2)+x2*x4p2*(x4-x2)+x1p2*(2.0*x2p2-x4p2));
			temp9 = (x2-x3)*(x2-x4)*(x3-x4)*x4;
			temp10 = x1p2*(2.0*x2p2-2.0*x2*x3+x4*(x3-x4));
			temp11 = x1*(-2.0*x2p2*x3-x3p2*x4+x4p3+x2*(2.0*x3p2+x3*x4-x4p2));
			A3 = (temp1+temp2+temp3+temp4+b*(temp5+temp6+temp7+temp8))/(12.0*x3*(temp9+temp10+temp11));
			
			temp1 = -12.0*a*x1*x2*(x1-x3)*(x2-x3)-3.0*ap4*(x1+x2-x3);
			temp2 = -6.0*x3*ap2*(x1p2+x1*(x2-x3)+x2*(x2-x3));
			temp3 = 4.0*ap3*(x1p2+x1*x2+x2p2-x3p2);
			temp4 = 12.0*x1*x2*(x1-x3)*(x2-x3)+3.0*bp3*(x1+x2-x3);
			temp5 = 6.0*b*x3*(x1p2+x1*x2+x2p2-x1*x3-x2*x3);
			temp6 = -4.0*bp2*(x1p2+x1*x2+x2p2-x3p2);
			temp7 = (x2-x3)*(x2-x4)*(x3-x4)*x4;
			temp8 = x1p2*(2.0*x2p2-2.0*x2*x3+(x3-x4)*x4);
			temp9 = x1*(-2.0*x2p2*x3-x3p2*x4+x4p3+x2*(2.0*x3p2+x3*x4-x4p2));
			A4 = (temp1+temp2+temp3+b*(temp4+temp5+temp6))/(12.0*(temp7+temp8+temp9));
			
			dom1 = A1*ddom1intdh(x1,x,fM,fS)+A2*ddom1intdh(x2,x,fM,fS)+A3*ddom1intdh(x3,x,fM,fS)+A4*ddom1intdh(x4,x,fM,fS);

			double dom2 = 0.0;
			ConHyperGeom fCHGn1(-0.25,1.5), fCHG1(0.25,1.5), fCHG3(0.75,2.5);
			GenHyperGeom fGHG1(0.5,1.0,1.75,2.25), fGHG3(1.5,2.0,2.75,3.25);
			
			temp1 = 105.0*pow(fS,1.5)/sqrt(fM-x);
			temp2 = 12.0*pow(fS,2.0)*0.919063*fCHGn1.Function(-2.0*xval);
			temp3 = -9.0*sqrt(2.0)*fS*0.906402*fCHG1.Function(-2.0*xval);
			temp4 = 2.0*(x-fM)*0.919063*fCHG3.Function(-2.0*xval);
			temp5 = 630.0*pow(2.0,0.25)*(x-fM)*pow(fS,2.0)*fGHG1.Function(-2.0*xval);
			temp6 = -32.0*pow(2.0,0.25)*pow((x-fM),3.0)*fGHG3.Function(-2.0*xval);
			temp7 = 630.0*pow(2.0,0.75)*sqrt(PI)*pow(fS,3.0);
			dom2 = (-sqrt(fM-x)*(temp1*(temp2+(x-fM)*(temp3+temp4))+temp5+temp6))/temp7;

			yval = dom1+dom2;
		}
	}
		return yval;
}

/*
* Returning values in groups - derived classes should define
* their own non-virtual function called within this functon
* which maps in to out w/o requiring a virtual function call
* everytime. Default behavior is just to map the virtual functions
* above.
*/
dArrayT& GreenwoodWilliamson::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl   = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		cout << "*** ERROR! Greenwood and Williamson potential unavailable.";
		*pddU++ = -1.0e10;
	}
	return(out);
}

dArrayT& GreenwoodWilliamson::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double xval, diff, yval;
	
		if (fS==0)
		{
			cout << "\n*** Bad SIGMA value in GreenwoodWilliamson.cpp.\n";
			throw eBadInputValue;
		}
		else
		{
			diff = *pl++ - fM;
			xval = pow(diff/(2.0*fS),2.0);
		
			if (diff > 0)
			{
				ModBessel k1(0.25), k3(0.75), k5(1.25);
				double f0 = sqrt(fS*(*pl++-fM))*exp(-xval);
				double f1 = 2.0*(pow(*pl++-fM,2.0) + 2.0*pow(fS,2.0))*k1.Function(xval);
				double f2 = -pow(*pl++-fM,2.0)*(k3.Function(xval) + k5.Function(xval));
		
				yval = f0*(f1 + f2);
			}
			else if (diff == 0)
				yval = 6.09752*pow(fS,3.0); 
			else
			{	/* decompose "negative domain": [h,mu] + [mu,inf]	*/
				/*													*/
				/* no closed form solution; use Gaussian quadrature */
				/* on the function:									*/
				/*			(s-h)^(1.5) phi(s)						*/
				double dom1=0.0;
				double a = *pl++;
				double b = fM;
				double x1, x2, x3, x4;
				double A1, A2, A3, A4;
				double temp1, temp2, temp3, temp4, temp5, temp6;
				double temp7, temp8, temp9, temp10, temp11;
				double ap2, ap3, ap4;
				double bp2, bp3;
				double x1p2, x1p3;
				double x2p2, x2p3;
				double x3p2;
				double x4p2, x4p3;
			
				x1 = (35.0*a-sqrt(525.0+70.0*sqrt(30.0))*(a-b)+35.0*b)/70.0;
				x2 = (35.0*a+sqrt(525.0+70.0*sqrt(30.0))*(a-b)+35.0*b)/70.0;
				x3 = (35.0*a-sqrt(525.0-70.0*sqrt(30.0))*(a-b)+35.0*b)/70.0;
				x4 = (35.0*a+sqrt(525.0-70.0*sqrt(30.0))*(a-b)+35.0*b)/70.0;
			
				ap2 = pow(a,2.0);
				ap3 = pow(a,3.0);
				ap4 = pow(a,4.0);
				bp2 = pow(b,2.0);
				bp3 = pow(b,3.0);
				x1p2 = pow(x1,2.0);
				x1p3 = pow(x1,3.0);
				x2p2 = pow(x2,2.0);
				x2p3 = pow(x2,3.0);
				x3p2 = pow(x3,2.0);
				x4p2 = pow(x4,2.0);
				x4p3 = pow(x4,3.0);
			
				temp1 = 12.0*a*x2*(x2-x3)*(x2-x4)*(x3-x4)*x4;
				temp2 = 4.0*ap3*(2.0*x2p3-2.0*x2*x3p2+x4*x3p2-x4p3);
				temp3 = ap4*(-6.0*x2p2+6.0*x2*x3+3.0*x4*(x4-x3));
				temp4 = 6.0*ap2*x3*(-2.0*x2p3+2.0*x2p2*x3+x4p2*(x4-x3));
				temp5 = -12.0*x2*(x2-x3)*(x2-x4)*(x3-x4)*x4;
				temp6 = 3*bp3*(2.0*x2p2-2.0*x2*x3+(x3-x4)*x4);
				temp7 = 6.0*b*x3*(2.0*x2p3-2.0*x2p2*x3+(x3-x4)*x4p2);
				temp8 = bp2*(-8.0*x2p3+8.0*x2*x3p2-4.0*x3p2*x4+4.0*x4p3);
				temp9 = (x2-x3)*(x2-x4)*(x3-x4)*x4;
				temp10 = x1p2*(2.0*x2p2-2.0*x2*x3+(x3-x4)*x4);
				temp11 = x1*(-2.0*x2p2*x3-x3p2*x4+x4p3+x2*(2.0*x3p2+x3*x4-x4p2));
				A1 = (temp1+temp2+temp3+temp4+b*(temp5+temp6+temp7+temp8))/(12.0*(x1-x2)*(temp9+temp10+temp11));
				
				temp1 = -12.0*a*x1*(x1-x3)*(x1-x4)*(x3-x4)*x4;
				temp2 = 3.0*ap4*(2.0*x1p2-2.0*x1*x3+(x3-x4)*x4);
				temp3 = 6.0*ap2*x3*(2.0*x1p3-2.0*x1p2*x3+(x3-x4)*x4p2);
				temp4 = ap3*(-8.0*x1p3+8.0*x1*x3p2-4.0*x3p2*x4+4.0*x4p3);
				temp5 = 12.0*x1*(x1-x3)*(x1-x4)*(x3-x4)*x4;
				temp6 = 4*bp2*(2.0*x1p3-2.0*x1*x3p2+x3p2*x4-x4p3);
				temp7 = bp3*(-6.0*x1p2+6.0*x1*x3+3.0*x4*(x4-x3));
				temp8 = 6.0*b*x3*(-2.0*x1p3+2.0*x1p2*x3+x4p2*(x4-x3));
				temp9 = (x2-x3)*(x2-x4)*(x3-x4)*x4;
				temp10 = x1p2*(2.0*x2p2-2.0*x2*x3+(x3-x4)*x4);
				temp11 = x1*(-2.0*x2p2*x3-x3p2*x4+x4p3+x2*(2.0*x3p2+x3*x4-x4p2));
				A2 = (temp1+temp2+temp3+temp4+b*(temp5+temp6+temp7+temp8))/(12.0*(x1-x2)*(temp9+temp10+temp11));
			
				temp1 = 12.0*a*x1*x2*(x1-x4)*(x2-x4)*x4;
				temp2 = -3.0*ap4*(2.0*x1*x2-x1*x4-x2*x4+x4p2);
				temp3 = 4.0*ap3*(x1p2*(2.0*x2-x4)+x1*x2*(2.0*x2-x4)-x2p2*x4+x4p3);
				temp4 = -6*ap2*(x1*x4p2*(x4-x2)+x2*x4p2*(x4-x2)+x1p2*(2.0*x2p2-x4p2));
				temp5 = -12.0*x1*x2*(x1-x4)*(x2-x4)*x4;
				temp6 = 3.0*bp3*(2.0*x1*x2-x1*x4-x2*x4+x4p2);
				temp7 = -4.0*bp2*(x1p2*(2.0*x2-x4)+x1*x2*(2.0*x2-x4)-x2p2*x4+x4p3);
				temp8 = 6.0*b*(x1*x4p2*(x4-x2)+x2*x4p2*(x4-x2)+x1p2*(2.0*x2p2-x4p2));
				temp9 = (x2-x3)*(x2-x4)*(x3-x4)*x4;
				temp10 = x1p2*(2.0*x2p2-2.0*x2*x3+x4*(x3-x4));
				temp11 = x1*(-2.0*x2p2*x3-x3p2*x4+x4p3+x2*(2.0*x3p2+x3*x4-x4p2));
				A3 = (temp1+temp2+temp3+temp4+b*(temp5+temp6+temp7+temp8))/(12.0*x3*(temp9+temp10+temp11));
			
				temp1 = -12.0*a*x1*x2*(x1-x3)*(x2-x3)-3.0*ap4*(x1+x2-x3);
				temp2 = -6.0*x3*ap2*(x1p2+x1*(x2-x3)+x2*(x2-x3));
				temp3 = 4.0*ap3*(x1p2+x1*x2+x2p2-x3p2);
				temp4 = 12.0*x1*x2*(x1-x3)*(x2-x3)+3.0*bp3*(x1+x2-x3);
				temp5 = 6.0*b*x3*(x1p2+x1*x2+x2p2-x1*x3-x2*x3);
				temp6 = -4.0*bp2*(x1p2+x1*x2+x2p2-x3p2);
				temp7 = (x2-x3)*(x2-x4)*(x3-x4)*x4;
				temp8 = x1p2*(2.0*x2p2-2.0*x2*x3+(x3-x4)*x4);
				temp9 = x1*(-2.0*x2p2*x3-x3p2*x4+x4p3+x2*(2.0*x3p2+x3*x4-x4p2));
				A4 = (temp1+temp2+temp3+b*(temp4+temp5+temp6))/(12.0*(temp7+temp8+temp9));
			
				dom1 = A1*dom1int(x1,*pl++,fM,fS)+A2*dom1int(x2,*pl++,fM,fS)+A3*dom1int(x3,*pl++,fM,fS)+A4*dom1int(x4,*pl++,fM,fS);

				double dom2 = 0.0;
				ConHyperGeom fCHG3(-0.75,0.5), fCHG1(-0.25,1.5);
				GenHyperGeom fGHG(0.5,1.0,1.75,2.25);
			
				temp1 = sqrt(fM-*pl++)/(5.0*fS*sqrt(2.0*PI));
				temp2 = 5.0*pow(2.0,0.25)*pow(fS,1.5)/sqrt(fM-*pl++);
				temp3 = fS*0.906402*fCHG3.Function(-2.0*xval);
				temp4 = sqrt(2.0)*(fM-*pl++)*0.919063*fCHG1.Function(-2.0*xval);
				temp5 = -2.0*pow((fM-*pl++),2.0)*fGHG.Function(-2.0*xval);
				dom2 = temp1*(temp2*(temp3+temp4)+temp5);

				yval = dom1+dom2;
			}
		}
		
		*pU++ = yval;
	}
	return(out);
}

dArrayT& GreenwoodWilliamson::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl  = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double xval, diff, yval;
	
		if (fS==0)
		{
			cout << "\n*** Bad SIGMA value in GreenwoodWilliamson.cpp.\n";
			throw eBadInputValue;
		}
		else
		{
			diff = *pl++ - fM;
			xval = pow(diff/(2.0*fS),2.0);
		
			if (diff > 0)
			{
				ModBessel k1(0.25), k3(0.75), k5(1.25), kn1(-0.25), k7(1.75), kn3(-0.75), k9(2.25);
				double f0 = sqrt(fS*(*pl++-fM))*exp(-xval);
				double f1 = 2.0*(pow(*pl++-fM,2.0) + 2.0*pow(fS,2.0))*k1.Function(xval);
				double f2 = -pow(*pl++-fM,2.0)*(k3.Function(xval) + k5.Function(xval));
			
				double df0 = -exp(-xval)*(pow(*pl++-fM,2.0)-pow(fS,2.0))/(2.0*pow(fS,1.5)*sqrt(*pl++-fM));
				double df1 = (*pl++-fM)*(4.0*k1.Function(xval)-(2.0*xval+1.0)*(kn3.Function(xval)+k5.Function(xval)));
				double df2sum = kn1.Function(xval)+k1.Function(xval)+k7.Function(xval)+k9.Function(xval);
				double df2 = (*pl++-fM)*(-2.0*(k3.Function(xval)+k5.Function(xval))+xval*df2sum);
	
				yval = df0*(f1+f2)+f0*(df1+df2);
			}
			else if (diff == 0)
				yval = -3.0*1.22542*sqrt(fS/PI)/pow(2.0,1.75);
			else
			{	/* decompose "negative domain": [h,mu] + [mu,inf]	*/
				/*													*/
				/* no closed form solution; use Gaussian quadrature */
				/* on the function:									*/
				/*			(s-h)^(1.5) phi(s)						*/
				double dom1=0.0;
				double a = *pl++;
				double b = fM;
				double x1, x2, x3, x4;
				double A1, A2, A3, A4;
				double temp1, temp2, temp3, temp4, temp5, temp6;
				double temp7, temp8, temp9, temp10, temp11;
				double ap2, ap3, ap4;
				double bp2, bp3;
				double x1p2, x1p3;
				double x2p2, x2p3;
				double x3p2;
				double x4p2, x4p3;
				
				x1 = (35.0*a-sqrt(525.0+70.0*sqrt(30.0))*(a-b)+35.0*b)/70.0;
				x2 = (35.0*a+sqrt(525.0+70.0*sqrt(30.0))*(a-b)+35.0*b)/70.0;
				x3 = (35.0*a-sqrt(525.0-70.0*sqrt(30.0))*(a-b)+35.0*b)/70.0;
				x4 = (35.0*a+sqrt(525.0-70.0*sqrt(30.0))*(a-b)+35.0*b)/70.0;
				
				ap2 = pow(a,2.0);
				ap3 = pow(a,3.0);
				ap4 = pow(a,4.0);
				bp2 = pow(b,2.0);
				bp3 = pow(b,3.0);
				x1p2 = pow(x1,2.0);
				x1p3 = pow(x1,3.0);
				x2p2 = pow(x2,2.0);
				x2p3 = pow(x2,3.0);
				x3p2 = pow(x3,2.0);
				x4p2 = pow(x4,2.0);
				x4p3 = pow(x4,3.0);
			
				temp1 = 12.0*a*x2*(x2-x3)*(x2-x4)*(x3-x4)*x4;
				temp2 = 4.0*ap3*(2.0*x2p3-2.0*x2*x3p2+x4*x3p2-x4p3);
				temp3 = ap4*(-6.0*x2p2+6.0*x2*x3+3.0*x4*(x4-x3));
				temp4 = 6.0*ap2*x3*(-2.0*x2p3+2.0*x2p2*x3+x4p2*(x4-x3));
				temp5 = -12.0*x2*(x2-x3)*(x2-x4)*(x3-x4)*x4;
				temp6 = 3*bp3*(2.0*x2p2-2.0*x2*x3+(x3-x4)*x4);
				temp7 = 6.0*b*x3*(2.0*x2p3-2.0*x2p2*x3+(x3-x4)*x4p2);
				temp8 = bp2*(-8.0*x2p3+8.0*x2*x3p2-4.0*x3p2*x4+4.0*x4p3);
				temp9 = (x2-x3)*(x2-x4)*(x3-x4)*x4;
				temp10 = x1p2*(2.0*x2p2-2.0*x2*x3+(x3-x4)*x4);
				temp11 = x1*(-2.0*x2p2*x3-x3p2*x4+x4p3+x2*(2.0*x3p2+x3*x4-x4p2));
				A1 = (temp1+temp2+temp3+temp4+b*(temp5+temp6+temp7+temp8))/(12.0*(x1-x2)*(temp9+temp10+temp11));
			
				temp1 = -12.0*a*x1*(x1-x3)*(x1-x4)*(x3-x4)*x4;
				temp2 = 3.0*ap4*(2.0*x1p2-2.0*x1*x3+(x3-x4)*x4);
				temp3 = 6.0*ap2*x3*(2.0*x1p3-2.0*x1p2*x3+(x3-x4)*x4p2);
				temp4 = ap3*(-8.0*x1p3+8.0*x1*x3p2-4.0*x3p2*x4+4.0*x4p3);
				temp5 = 12.0*x1*(x1-x3)*(x1-x4)*(x3-x4)*x4;
				temp6 = 4*bp2*(2.0*x1p3-2.0*x1*x3p2+x3p2*x4-x4p3);
				temp7 = bp3*(-6.0*x1p2+6.0*x1*x3+3.0*x4*(x4-x3));
				temp8 = 6.0*b*x3*(-2.0*x1p3+2.0*x1p2*x3+x4p2*(x4-x3));
				temp9 = (x2-x3)*(x2-x4)*(x3-x4)*x4;
				temp10 = x1p2*(2.0*x2p2-2.0*x2*x3+(x3-x4)*x4);
				temp11 = x1*(-2.0*x2p2*x3-x3p2*x4+x4p3+x2*(2.0*x3p2+x3*x4-x4p2));
				A2 = (temp1+temp2+temp3+temp4+b*(temp5+temp6+temp7+temp8))/(12.0*(x1-x2)*(temp9+temp10+temp11));
			
				temp1 = 12.0*a*x1*x2*(x1-x4)*(x2-x4)*x4;
				temp2 = -3.0*ap4*(2.0*x1*x2-x1*x4-x2*x4+x4p2);
				temp3 = 4.0*ap3*(x1p2*(2.0*x2-x4)+x1*x2*(2.0*x2-x4)-x2p2*x4+x4p3);
				temp4 = -6*ap2*(x1*x4p2*(x4-x2)+x2*x4p2*(x4-x2)+x1p2*(2.0*x2p2-x4p2));
				temp5 = -12.0*x1*x2*(x1-x4)*(x2-x4)*x4;
				temp6 = 3.0*bp3*(2.0*x1*x2-x1*x4-x2*x4+x4p2);
				temp7 = -4.0*bp2*(x1p2*(2.0*x2-x4)+x1*x2*(2.0*x2-x4)-x2p2*x4+x4p3);
				temp8 = 6.0*b*(x1*x4p2*(x4-x2)+x2*x4p2*(x4-x2)+x1p2*(2.0*x2p2-x4p2));
				temp9 = (x2-x3)*(x2-x4)*(x3-x4)*x4;
				temp10 = x1p2*(2.0*x2p2-2.0*x2*x3+x4*(x3-x4));
				temp11 = x1*(-2.0*x2p2*x3-x3p2*x4+x4p3+x2*(2.0*x3p2+x3*x4-x4p2));
				A3 = (temp1+temp2+temp3+temp4+b*(temp5+temp6+temp7+temp8))/(12.0*x3*(temp9+temp10+temp11));
			
				temp1 = -12.0*a*x1*x2*(x1-x3)*(x2-x3)-3.0*ap4*(x1+x2-x3);
				temp2 = -6.0*x3*ap2*(x1p2+x1*(x2-x3)+x2*(x2-x3));
				temp3 = 4.0*ap3*(x1p2+x1*x2+x2p2-x3p2);
				temp4 = 12.0*x1*x2*(x1-x3)*(x2-x3)+3.0*bp3*(x1+x2-x3);
				temp5 = 6.0*b*x3*(x1p2+x1*x2+x2p2-x1*x3-x2*x3);
				temp6 = -4.0*bp2*(x1p2+x1*x2+x2p2-x3p2);
				temp7 = (x2-x3)*(x2-x4)*(x3-x4)*x4;
				temp8 = x1p2*(2.0*x2p2-2.0*x2*x3+(x3-x4)*x4);
				temp9 = x1*(-2.0*x2p2*x3-x3p2*x4+x4p3+x2*(2.0*x3p2+x3*x4-x4p2));
				A4 = (temp1+temp2+temp3+b*(temp4+temp5+temp6))/(12.0*(temp7+temp8+temp9));
			
				dom1 = A1*ddom1intdh(x1,*pl++,fM,fS)+A2*ddom1intdh(x2,*pl++,fM,fS)+A3*ddom1intdh(x3,*pl++,fM,fS)+A4*ddom1intdh(x4,*pl++,fM,fS);
	
				double dom2 = 0.0;
				ConHyperGeom fCHGn1(-0.25,1.5), fCHG1(0.25,1.5), fCHG3(0.75,2.5);
				GenHyperGeom fGHG1(0.5,1.0,1.75,2.25), fGHG3(1.5,2.0,2.75,3.25);
			
				temp1 = 105.0*pow(fS,1.5)/sqrt(fM-*pl++);
				temp2 = 12.0*pow(fS,2.0)*0.919063*fCHGn1.Function(-2.0*xval);
				temp3 = -9.0*sqrt(2.0)*fS*0.906402*fCHG1.Function(-2.0*xval);
				temp4 = 2.0*(*pl++-fM)*0.919063*fCHG3.Function(-2.0*xval);
				temp5 = 630.0*pow(2.0,0.25)*(*pl++-fM)*pow(fS,2.0)*fGHG1.Function(-2.0*xval);
				temp6 = -32.0*pow(2.0,0.25)*pow((*pl++-fM),3.0)*fGHG3.Function(-2.0*xval);
				temp7 = 630.0*pow(2.0,0.75)*sqrt(PI)*pow(fS,3.0);
				dom2 = (-sqrt(fM-*pl++)*(temp1*(temp2+(*pl++-fM)*(temp3+temp4))+temp5+temp6))/temp7;

				yval = dom1+dom2;
			}
		}
	
		*pdU++ = yval;
	}
	return(out);
}





