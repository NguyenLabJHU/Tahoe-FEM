/* $Id: CubicSplineT.cpp,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (12/02/1996)                                          */
/* CubicSplineT.cpp                                                       */

#include "CubicSplineT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"
#include "TriDiagdMatrixT.h"
#include "iArrayT.h"

const int kNumSplineCoeffs = 4;

/* constructor */
CubicSplineT::CubicSplineT(const dArrayT& knots, const dArray2DT& coefficients):
	fXPoints(knots),
	fCoefficients(coefficients)
{
	/* check dimensions */
	if (fCoefficients.MajorDim() != fXPoints.Length() + 1)
	    throw eSizeMismatch;
	if (fCoefficients.MinorDim() != kNumSplineCoeffs)
	    throw eGeneralFail;
}

CubicSplineT::CubicSplineT(const dArray2DT& points, FixityT fixity):
	fXPoints(0,points),
	fCoefficients(fXPoints.Length() + 1, kNumSplineCoeffs)
{
	/* compute spline coefficients */
	SetSpline(points, fixity);
}

/* I/O */
void CubicSplineT::Print(ostream& out) const
{
	/* parameters */
	out << " Knots:\n\n";
	out << fXPoints;
}

void CubicSplineT::PrintName(ostream& out) const
{
	out << "    Cubic spline\n";
}
	    	
/* returning values */
double CubicSplineT::Function(double x) const { return function(x); }
double CubicSplineT::DFunction(double x) const { return Dfunction(x); }
double CubicSplineT::DDFunction(double x) const { return DDfunction(x); }

/* returning values in groups - returns refence to out to allow:
*
*	dArrayT& goodname = pfunc->MapFunction(in, tempspace);
*/
dArrayT& CubicSplineT::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension check */
	if ( in.Length() != out.Length() ) throw eGeneralFail;
	
	double *pin   = in.Pointer();
	double *pout  = out.Pointer();
	int    length = in.Length();
	
	/* fast mapping */
	for (int i = 0; i < length; i++)
		*pout++ = function(*pin++);

	return out;
}

dArrayT& CubicSplineT::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension check */
	if ( in.Length() != out.Length() ) throw eGeneralFail;
	
	double *pin   =  in.Pointer();
	double *pout  = out.Pointer();
	int    length = in.Length();
	
	/* fast mapping */
	for (int i = 0; i < length; i++)
		*pout++ = Dfunction(*pin++);

	return out;
}

dArrayT& CubicSplineT::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension check */
	if ( in.Length() != out.Length() ) throw eGeneralFail;
	
	double *pin   =  in.Pointer();
	double *pout  = out.Pointer();
	int    length = in.Length();
	
	/* fast mapping */
	for (int i = 0; i < length; i++)
		*pout++ = DDfunction(*pin++);

	return out;
}

/*
* Return 0th, 1st, and 2nd derivative in the respective
* fields of the dArrayT.
*/  	
void CubicSplineT::SetAll(double x, dArrayT& data) const
{
	all_functions(x, data[0], data[1], data[2]);
}

/* accessor to the spline information */
const dArray2DT& CubicSplineT::Coefficients(void) const
{
	return fCoefficients;
}

/**********************************************************************
* Protected
**********************************************************************/

/* non-virtual function calls */
double CubicSplineT::function(double x) const
{
	int i = fXPoints.Range(x);
	
	double  dx = (i == 0) ? x - fXPoints[i] :
	             ((i == fXPoints.Length()) ? x - fXPoints[i-1] :
	              x - fXPoints[i-1]);
	double*  a = fCoefficients(i);	
	
	return a[0] + a[1]*dx + a[2]*dx*dx + a[3]*dx*dx*dx;
}

double CubicSplineT::Dfunction(double x) const
{
	int i = fXPoints.Range(x);
	
	double  dx = (i == 0) ? x - fXPoints[i] :
	             ((i == fXPoints.Length()) ? x - fXPoints[i-1] :
	              x - fXPoints[i-1]);
	double*  a = fCoefficients(i);	
	
	return a[1] + 2.0*a[2]*dx + 3.0*a[3]*dx*dx;
}

double CubicSplineT::DDfunction(double x) const	
{
	int i = fXPoints.Range(x);
	
	double  dx = (i == 0) ? x - fXPoints[i] :
	             ((i == fXPoints.Length()) ? x - fXPoints[i-1] :
	              x - fXPoints[i-1]);
	double*  a = fCoefficients(i);	
	
	return 2.0*a[2] + 6.0*a[3]*dx;
}

void CubicSplineT::all_functions(double x, double& f, double& Df, double& DDf) const
{
	int i = fXPoints.Range(x);
	
	double  dx = (i == 0) ? x - fXPoints[i] :
	             ((i == fXPoints.Length()) ? x - fXPoints[i-1] :
	              x - fXPoints[i-1]);
	double*  a = fCoefficients(i);	
	
	f = a[0] + a[1]*dx + a[2]*dx*dx + a[3]*dx*dx*dx;
	Df = a[1] + 2.0*a[2]*dx + 3.0*a[3]*dx*dx;
	DDf = 2.0*a[2] + 6.0*a[3]*dx;
}

/**********************************************************************
* Private
**********************************************************************/

/* compute spline coefficients */
void CubicSplineT::SetSpline(const dArray2DT& points, FixityT fixity)
{
	/* allocate work space */
	dArrayT YPoints(fXPoints.Length());
	dArrayT DDY(fXPoints.Length());
	dArrayT dxi(fXPoints.Length() - 1);
	
	/* copy y-values */
	points.ColumnCopy(1,YPoints);
	
	/* set delta's */
	for (int ii = 1; ii < fXPoints.Length(); ii++)
		dxi[ii-1] = fXPoints[ii] - fXPoints[ii-1];
	
	/* construct solver */
	int numeqs = fXPoints.Length() - 2;
	TriDiagdMatrixT LHS(numeqs);
	dArrayT RHS(numeqs, DDY.Pointer() + 1);

	/* build equation system */
	double* pdx  = dxi.Pointer();
	double* pY   = YPoints.Pointer() + 1;
	double* pRHS = RHS.Pointer();
	for (int i = 0; i < numeqs; i++)
	{
		LHS.SetRow(i, *pdx/6.0, (*pdx + *(pdx+1))/3.0, *(pdx+1)/6.0);
		*pRHS = ((*(pY+1) - *pY)/(*(pdx+1))) -
		        ((*pY - *(pY-1))/(*pdx));
		
		pdx++; pY++; pRHS++;
	}	

	/* set end conditions */
	if (fixity == kParabolic)
	{
		LHS.AddToRow(0,0.0,dxi[0]/6.0,0.0);
		LHS.AddToRow(numeqs-1,0.0,dxi[numeqs-1]/6.0,0.0);
	}
	
	/* solve */
	LHS.LinearSolve(RHS);

	/* end conditions */
	if (fixity == kFreeRun)
		DDY[0] = DDY[numeqs+1] = 0.0;
	else if (fixity == kParabolic)
	{
		DDY[0] = DDY[1];
		DDY[numeqs+1] = DDY[numeqs];
	}
	else
		throw eBadInputValue;
	
	/* store spline coefficients */
	for (int j = 1; j < fXPoints.Length(); j++)
	{
		int i = j-1;
		double* coeff = fCoefficients(j);
		double dx = dxi[i];
		
		coeff[0] = YPoints[i];
		coeff[1] = -dx*(2.0*DDY[i] + DDY[i+1])/6.0 +
		               (YPoints[i+1] - YPoints[i])/dx;
		coeff[2] = DDY[i]/2.0;
		coeff[3] = (DDY[i+1] - DDY[i])/(6.0*dx);
	}
	
	/* extensions */
	fCoefficients(0,0) = fCoefficients(1,0);
	fCoefficients(0,1) = fCoefficients(1,1);
	fCoefficients(0,2) = fCoefficients(1,2);
	fCoefficients(0,3) = 0.0;

	int dex = fXPoints.Length() - 1;
	fCoefficients(dex+1,0) = YPoints[dex];
	fCoefficients(dex+1,1) = dxi[dex-1]*(DDY[dex-1] + 2.0*DDY[dex])/6.0 +
		               (YPoints[dex] - YPoints[dex-1])/dxi[dex-1];
	fCoefficients(dex+1,2) = DDY[dex]/2.0;
	fCoefficients(dex+1,3) = 0.0;
}
