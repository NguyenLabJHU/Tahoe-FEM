/* $Id: CosinePlusT.cpp,v 1.1 2009-05-20 17:48:26 regueiro Exp $ */
#include "CosinePlusT.h"
#include "dArrayT.h"
#include <math.h>

using namespace Tahoe;

/* constructor */
CosinePlusT::CosinePlusT(double a, double b, double c, double d, double e, double f, double g,
		double p, double q):
	fa(a),
	fb(b),
	fc(c),
	fd(d),
	fe(e),
	ff(f),
	fg(g),
	fp(p),
	fq(q)
{
	SetName("cosine_plus");
}

CosinePlusT::CosinePlusT(void):
	fa(0.0),
	fb(0.0),
	fc(0.0),
	fd(0.0),
	fe(0.0),
	ff(0.0),
	fg(0.0),
	fp(0.0),
	fq(0.0)
{
	SetName("cosine_plus");
}

/** evaluate function */
double CosinePlusT::Function(double x) const
{
	return fa + fb*cos(fc*x) + fd*sin(fe*x) + ff*x*cos(fg*x) + fp*x*cos(fq*x);
}

/** evaluate first derivative function */
double CosinePlusT::DFunction(double x) const
{
	return -fb*fc*sin(fc*x) + fd*fe*cos(fe*x) + ff*cos(fg*x) - ff*fg*x*sin(fg*x) + fp*sin(fq*x) + fp*fq*x*cos(fq*x);
}

/** evaluate second derivative function */
double CosinePlusT::DDFunction(double x) const
{
	return -fb*fc*fc*cos(fc*x) - fd*fe*fe*sin(fe*x) - 2*ff*fg*sin(fg*x) - ff*fg*fg*x*cos(fg*x) + 2*fp*fq*cos(fq*x) - fp*fq*fq*x*sin(fq*x);
}

/* Returning values in groups */

/** multiple function evaluations */
dArrayT& CosinePlusT::MapFunction(const dArrayT& in, dArrayT& out) const
{
#if __option(extended_errorcheck)
	/* dimension checks */
	if (in.Length() != out.Length()) ExceptionT::GeneralFail("CosinePlusT::MapFunction(");
#endif

	int length = in.Length();
	double* y = out.Pointer();
	const double* x = in.Pointer();
	for (int i = 0; i < length; i++)
	{
		double r = *x++;
		
		*y++ = fa + fb*cos(fc*r) + fd*sin(fe*r) + ff*r*cos(fg*r) + fp*r*cos(fq*r);
	}

	return out;
}

/** multiple first derivative evaluations */
dArrayT& CosinePlusT::MapDFunction(const dArrayT& in, dArrayT& out) const
{
#if __option(extended_errorcheck)
	/* dimension checks */
	if (in.Length() != out.Length()) ExceptionT::GeneralFail("CosinePlusT::MapDFunction(");
#endif

	int length = in.Length();
	double* y = out.Pointer();
	const double* x = in.Pointer();
	for (int i = 0; i < length; i++)
	{
		double r = *x++;
		
		*y++ = -fb*fc*sin(fc*r) + fd*fe*cos(fe*r) + ff*cos(fg*r) - ff*fg*r*sin(fg*r) + fp*sin(fq*r) + fp*fq*r*cos(fq*r);
	}

	return out;
}

/** multiple second derivative evaluations */
dArrayT& CosinePlusT::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
#if __option(extended_errorcheck)
	/* dimension checks */
	if (in.Length() != out.Length()) ExceptionT::GeneralFail("CosinePlusT::MapDDFunction(");
#endif

	int length = in.Length();
	double* y = out.Pointer();
	const double* x = in.Pointer();
	for (int i = 0; i < length; i++)
	{
		double r = *x++;
		
		*y++ = -fb*fc*fc*cos(fc*r) - fd*fe*fe*sin(fe*r) - 2*ff*fg*sin(fg*r) - ff*fg*fg*r*cos(fg*r) + 2*fp*fq*cos(fq*r) - fp*fq*fq*r*sin(fq*r);
	}

	return out;
}

/* describe the parameters needed by the interface */
void CosinePlusT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);
	
	list.AddParameter(fa, "a");
	list.AddParameter(fb, "b");
	list.AddParameter(fc, "c");
	list.AddParameter(fd, "d");
	list.AddParameter(fe, "e");
	list.AddParameter(ff, "f");
	list.AddParameter(fg, "g");
	list.AddParameter(fp, "p");
	list.AddParameter(fq, "q");
	
	/* set the description */
	list.SetDescription("f(x) = a + b cos(c t) + d sin(e t) + f t cos(g t) + p t sin(q t)");
}

/* accept parameter list */
void CosinePlusT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);
	
	fa = list.GetParameter("a");
	fb = list.GetParameter("b");
	fc = list.GetParameter("c");
	fd = list.GetParameter("d");
	fe = list.GetParameter("e");
	ff = list.GetParameter("f");
	fg = list.GetParameter("g");
	fp = list.GetParameter("p");
	fq = list.GetParameter("q");
}
