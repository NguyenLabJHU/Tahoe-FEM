/* $Id: CosineT.cpp,v 1.1 2004-06-29 16:15:43 paklein Exp $ */
#include "CosineT.h"
#include "dArrayT.h"
#include <math.h>

using namespace Tahoe;

/* constructor */
CosineT::CosineT(double a, double b, double c, double d):
	fa(a),
	fb(b),
	fc(c),
	fd(d)
{
	SetName("cosine");
}

CosineT::CosineT(void):
	fa(0.0),
	fb(0.0),
	fc(0.0),
	fd(0.0)
{
	SetName("cosine");
}

/** evaluate function */
double CosineT::Function(double x) const
{
	return fa + fb*cos(fc*x + fd);
}

/** evaluate first derivative function */
double CosineT::DFunction(double x) const
{
	return -fb*fc*sin(fc*x + fd);
}

/** evaluate second derivative function */
double CosineT::DDFunction(double x) const
{
	return -fb*fc*fc*cos(fc*x + fd);
}

/* Returning values in groups */

/** multiple function evaluations */
dArrayT& CosineT::MapFunction(const dArrayT& in, dArrayT& out) const
{
#if __option(extended_errorcheck)
	/* dimension checks */
	if (in.Length() != out.Length()) ExceptionT::GeneralFail("CosineT::MapFunction(");
#endif

	int length = in.Length();
	double* y = out.Pointer();
	const double* x = in.Pointer();
	for (int i = 0; i < length; i++)
		*y++ = fa + fb*cos(fc*(*x++) + fd);

	return out;
}

/** multiple first derivative evaluations */
dArrayT& CosineT::MapDFunction(const dArrayT& in, dArrayT& out) const
{
#if __option(extended_errorcheck)
	/* dimension checks */
	if (in.Length() != out.Length()) ExceptionT::GeneralFail("CosineT::MapDFunction(");
#endif

	int length = in.Length();
	double* y = out.Pointer();
	const double* x = in.Pointer();
	for (int i = 0; i < length; i++)
		*y++ = -fb*fc*sin(fc*(*x++) + fd);

	return out;
}

/** multiple second derivative evaluations */
dArrayT& CosineT::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
#if __option(extended_errorcheck)
	/* dimension checks */
	if (in.Length() != out.Length()) ExceptionT::GeneralFail("CosineT::MapDDFunction(");
#endif

	int length = in.Length();
	double* y = out.Pointer();
	const double* x = in.Pointer();
	for (int i = 0; i < length; i++)
		*y++ = -fb*fc*fc*cos(fc*(*x++) + fd);

	return out;
}

/* describe the parameters needed by the interface */
void CosineT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);
	
	list.AddParameter(fa, "a");
	list.AddParameter(fb, "b");
	list.AddParameter(fc, "c");
	list.AddParameter(fd, "d");
	
	/* set the description */
	list.SetDescription("f(x) = a + b cos(c t + d)");
}

/* accept parameter list */
void CosineT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);
	
	fa = list.GetParameter("a");
	fb = list.GetParameter("b");
	fc = list.GetParameter("c");
	fd = list.GetParameter("d");
}
