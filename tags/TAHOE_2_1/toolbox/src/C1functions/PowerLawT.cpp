/* $Id: PowerLawT.cpp,v 1.1 2004-01-27 19:07:23 paklein Exp $ */
#include "PowerLawT.h"
#include "dArrayT.h"
#include <math.h>

using namespace Tahoe;

/* constructor */
PowerLawT::PowerLawT(double a, double b, double n):
	fa(a),
	fb(b),
	fn(n)
{
	SetName("power_law");
}

PowerLawT::PowerLawT(void):
	fa(0.0),
	fb(0.0),
	fn(0.0)
{
	SetName("power_law");
}

/** evaluate function */
double PowerLawT::Function(double x) const
{
	return fa*pow(1.0 + fb*x, fn);
}

/** evaluate first derivative function */
double PowerLawT::DFunction(double x) const
{
	return fa*fb*fn*pow(1.0 + fb*x, fn - 1.0);
}

/** evaluate second derivative function */
double PowerLawT::DDFunction(double x) const
{
	return fa*fb*fb*fn*(fn - 1.0)*pow(1.0 + fb*x, fn - 2.0);
}

/* Returning values in groups */

/** multiple function evaluations */
dArrayT& PowerLawT::MapFunction(const dArrayT& in, dArrayT& out) const
{
#if __option(extended_errorcheck)
	/* dimension checks */
	if (in.Length() != out.Length()) ExceptionT::GeneralFail("PowerLawT::MapFunction(");
#endif

	int length = in.Length();
	double* y = out.Pointer();
	const double* x = in.Pointer();
	for (int i = 0; i < length; i++)
		*y++ = fa*pow(1.0 + fb*(*x++), fn);

	return out;
}

/** multiple first derivative evaluations */
dArrayT& PowerLawT::MapDFunction(const dArrayT& in, dArrayT& out) const
{
#if __option(extended_errorcheck)
	/* dimension checks */
	if (in.Length() != out.Length()) ExceptionT::GeneralFail("PowerLawT::MapDFunction(");
#endif

	int length = in.Length();
	double* y = out.Pointer();
	const double* x = in.Pointer();
	for (int i = 0; i < length; i++)
		*y++ = fa*fb*fn*pow(1.0 + fb*(*x++), fn - 1.0);

	return out;
}

/** multiple second derivative evaluations */
dArrayT& PowerLawT::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
#if __option(extended_errorcheck)
	/* dimension checks */
	if (in.Length() != out.Length()) ExceptionT::GeneralFail("PowerLawT::MapDDFunction(");
#endif

	int length = in.Length();
	double* y = out.Pointer();
	const double* x = in.Pointer();
	for (int i = 0; i < length; i++)
		*y++ = fa*fb*fb*fn*(fn - 1.0)*pow(1.0 + fb*(*x++), fn - 2.0);

	return out;
}

/* describe the parameters needed by the interface */
void PowerLawT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);
	
	list.AddParameter(fa, "a");
	list.AddParameter(fb, "b");
	list.AddParameter(fn, "n");
	
	/* set the description */
	list.SetDescription("f(x) = a*(1 + b*x)^n");
}

/* accept parameter list */
void PowerLawT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);
	
	fa = list.GetParameter("a");
	fb = list.GetParameter("b");
	fn = list.GetParameter("n");
}
