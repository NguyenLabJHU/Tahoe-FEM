/* $Id: PiecewiseLinearT.cpp,v 1.4 2004-12-27 06:07:39 paklein Exp $ */
#include "PiecewiseLinearT.h"
#include "dArray2DT.h"

using namespace Tahoe;

/* constructor */
PiecewiseLinearT::PiecewiseLinearT(void)
{
	SetName("piecewise_linear");
}

PiecewiseLinearT::PiecewiseLinearT(const dArray2DT& points)
{
	SetName("piecewise_linear");
	SetPoints(points);
}

/* set function values */
void PiecewiseLinearT::SetPoints(const dArray2DT& points)
{
	/* collect data */
	fXPoints.SetValues(0, points);
	fYPoints.Dimension(fXPoints.Length());
	points.ColumnCopy(1, fYPoints);

	/* check */
	CheckKnots(fXPoints, fYPoints);
}

/* I/O */
void PiecewiseLinearT::Print(ostream& out) const
{
	int d_width = OutputWidth(out, &fXPoints[0]);

	/* parameters */
	out << " Points:\n";
	for (int i = 0; i < fXPoints.Length(); i++)
		out << setw(kIntWidth) << i+1 
		    << setw(d_width) << fXPoints[i] 
		    << setw(d_width) << Function(fXPoints[i]) << '\n';
}

void PiecewiseLinearT::PrintName(ostream& out) const
{
	out << "    Piecewise linear\n";
}
	    	
/* returning values */
double PiecewiseLinearT::Function(double x) const { return function(x); }
double PiecewiseLinearT::DFunction(double x) const { return Dfunction(x); }
double PiecewiseLinearT::DDFunction(double x) const { return DDfunction(x); }
double PiecewiseLinearT::DDDFunction(double x) const { return DDfunction(x); }
double PiecewiseLinearT::DDDDFunction(double x) const { return DDfunction(x); }

/* returning values in groups - returns refence to out to allow:
 *
 *	dArrayT& goodname = pfunc->MapFunction(in, tempspace);
 */
dArrayT& PiecewiseLinearT::MapFunction(const dArrayT& in, dArrayT& out) const
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (in.Length() != out.Length()) ExceptionT::SizeMismatch("PiecewiseLinearT::MapFunction");
#endif
	
	const double *pin   = in.Pointer();
	double *pout  = out.Pointer();
	int    length = in.Length();
	
	/* fast mapping */
	for (int i = 0; i < length; i++)
		*pout++ = function(*pin++);

	return out;
}

dArrayT& PiecewiseLinearT::MapDFunction(const dArrayT& in, dArrayT& out) const
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (in.Length() != out.Length()) ExceptionT::SizeMismatch("PiecewiseLinearT::MapDFunction");
#endif
	
	const double *pin   =  in.Pointer();
	double *pout  = out.Pointer();
	int    length = in.Length();
	
	/* fast mapping */
	for (int i = 0; i < length; i++)
		*pout++ = Dfunction(*pin++);

	return out;
}

dArrayT& PiecewiseLinearT::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (in.Length() != out.Length()) ExceptionT::SizeMismatch("PiecewiseLinearT::MapDDFunction");
#endif
	
	out = 0.0;
	return out;
}

dArrayT& PiecewiseLinearT::MapDDDFunction(const dArrayT& in, dArrayT& out) const
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (in.Length() != out.Length()) ExceptionT::SizeMismatch("PiecewiseLinearT::MapDDFunction");
#endif
	
	out = 0.0;
	return out;
}

dArrayT& PiecewiseLinearT::MapDDDDFunction(const dArrayT& in, dArrayT& out) const
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (in.Length() != out.Length()) ExceptionT::SizeMismatch("PiecewiseLinearT::MapDDFunction");
#endif
	
	out = 0.0;
	return out;
}

/*
* Return 0th, 1st, and 2nd derivative in the respective
* fields of the dArrayT.
*/  	
void PiecewiseLinearT::SetAll(double x, dArrayT& data) const
{
	all_functions(x, data[0], data[1], data[2]);
}

/* information about subordinate parameter lists */
void PiecewiseLinearT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	C1FunctionT::DefineSubs(sub_list);
	
	/* spline points */
	sub_list.AddSub("OrderedPair", ParameterListT::OnePlus);
}

/* accept parameter list */
void PiecewiseLinearT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);

	/* collect spline points */
	dArray2DT points(list.NumLists("OrderedPair"), 2);
	for (int i = 0; i < points.MajorDim(); i++)
	{
		const ParameterListT* knot = list.List("OrderedPair", i);
		points(i,0) = knot->GetParameter("x");
		points(i,1) = knot->GetParameter("y");
	}

	/* dimension internal data structures */
	dArrayT x_points(points.MajorDim());
	points.ColumnCopy(0, x_points);	
	fXPoints.SetValues(x_points);
	fYPoints.Dimension(fXPoints.Length());
	points.ColumnCopy(1, fYPoints);	

	/* check */
	CheckKnots(fXPoints, fYPoints);
}

/**********************************************************************
 * protected
 **********************************************************************/

/* non-virtual function calls */
double PiecewiseLinearT::function(double x) const
{
	int i = fXPoints.Range(x);
	if (i == 0) /* first abscissa */
		return fYPoints[0];
	else if (i == fXPoints.Length()) /* last abscissa */
		return fYPoints.Last();
	else /* linearly interpolate */
		return fYPoints[i-1] + (x - fXPoints[i-1])*(fYPoints[i] - fYPoints[i-1])/(fXPoints[i] - fXPoints[i-1]);
}

double PiecewiseLinearT::Dfunction(double x) const
{
	int i = fXPoints.Range(x);
	if (i == 0 || i == fXPoints.Length()) /* first of last abscissa */
		return 0.0;
	else /* slope over the interface */
		return (fYPoints[i] - fYPoints[i-1])/(fXPoints[i] - fXPoints[i-1]);
}

double PiecewiseLinearT::DDfunction(double x) const	
{ 
#pragma unused(x)
	return 0.0;
}

double PiecewiseLinearT::DDDfunction(double x) const	
{ 
#pragma unused(x)
	return 0.0;
}

double PiecewiseLinearT::DDDDfunction(double x) const	
{ 
#pragma unused(x)
	return 0.0;
}

void PiecewiseLinearT::all_functions(double x, double& f, double& Df, double& DDf) const
{
	int i = fXPoints.Range(x);
	double slope = 0.0;
	if (i > 0 && i < fXPoints.Length())
		slope = (fYPoints[i] - fYPoints[i-1])/(fXPoints[i] - fXPoints[i-1]);

	f   = fYPoints[i-1] + (x - fXPoints[i-1])*slope;
	Df  = slope;
	DDf = 0.0;
}

/**********************************************************************
 * Private
 **********************************************************************/

/* make sure knots are valid */
void PiecewiseLinearT::CheckKnots(const dRangeArrayT& x, const dArrayT& y) const
{
#pragma unused(y)
	for (int i = 1; i < x.Length(); i++)
		if (fabs(x[i] - x[i-1]) < kSmall)
			ExceptionT::GeneralFail("PiecewiseLinearT::CheckKnots", "knots %d and %d are coincident at %g",
				i, i+1, x[i]);
}
