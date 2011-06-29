/* $Id: LinearExponentialT.cpp,v 1.6 2004-07-20 23:23:33 rdorgan Exp $ */
/* created: paklein (10/30/1997) */
#include "LinearExponentialT.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructor */
LinearExponentialT::LinearExponentialT(double a, double b, double c, double d):
	fa(a),
	fb(b),
	fc(c),
	fd(d)
{
	SetName("linear_exponential");
	if (fabs(fd) < kSmall) throw ExceptionT::kBadInputValue;
}

LinearExponentialT::LinearExponentialT(void):
	fa(0.0), fb(0.0), fc(0.0), fd(0.0)
{
	SetName("linear_exponential");
}

/* print parameters */
void LinearExponentialT::Print(ostream& out) const
{
	/* parameters */
	out << " the function:\n";
	out << "     f(x) = a + b x + c (1 - exp[-d x])\n";
	out << " with:\n";
	out << "     a = " << fa << '\n';
	out << "     b = " << fb << '\n';
	out << "     c = " << fc << '\n';
	out << "     d = " << fd << '\n';
}

/* print function name */
void LinearExponentialT::PrintName(ostream& out) const
{
	out << "    Linear-exponential\n";
}

/* evaluate function */
double LinearExponentialT::Function(double x) const
{
	return fa + fb*x + fc*(1.0 - exp(-x/fd));
}

/* evaluate first derivative function */
double LinearExponentialT::DFunction(double x) const
{
	return fb + fc*exp(-x/fd)/fd;
}

/* evaluate second derivative function */
double LinearExponentialT::DDFunction(double x) const
{
	return -fc*exp(-x/fd)/fd/fd;
}

/* evaluate third derivative function */
double LinearExponentialT::DDDFunction(double x) const
{
	return fc*exp(-x/fd)/fd/fd/fd;
}

/* evaluate fourth derivative function */
double LinearExponentialT::DDDDFunction(double x) const
{
	return -fc*exp(-x/fd)/fd/fd/fd/fd;
}

/* Returning values in groups */

/* multiple function evaluations */
dArrayT& LinearExponentialT::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	int length = in.Length();
	double* y = out.Pointer();
	const double* x = in.Pointer();
	for (int i = 0; i < length; i++)
	{
		*y++ = fa + fb*(*x) + fc*(1.0 - exp(-(*x)/fd));
		x++;
	}
	return out;
}

/* multiple first derivative evaluations */
dArrayT& LinearExponentialT::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	int length = in.Length();
	double* y = out.Pointer();
	const double* x = in.Pointer();
	for (int i = 0; i < length; i++)
		*y++ = fb + fc*exp(-(*x++)/fd)/fd;
	return out;
}

/* multiple second derivative evaluations */
dArrayT& LinearExponentialT::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	int length = in.Length();
	double* y = out.Pointer();
	const double* x = in.Pointer();
	for (int i = 0; i < length; i++)
		*y++ = -fc*exp(-(*x++)/fd)/fd/fd;
	return out;
}

/* multiple third derivative evaluations */
dArrayT& LinearExponentialT::MapDDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	int length = in.Length();
	double* y = out.Pointer();
	const double* x = in.Pointer();
	for (int i = 0; i < length; i++)
		*y++ = fc*exp(-(*x++)/fd)/fd/fd/fd;
	return out;
}

/* multiple fourth derivative evaluations */
dArrayT& LinearExponentialT::MapDDDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	int length = in.Length();
	double* y = out.Pointer();
	const double* x = in.Pointer();
	for (int i = 0; i < length; i++)
		*y++ = -fc*exp(-(*x++)/fd)/fd/fd/fd/fd;
	return out;
}

void LinearExponentialT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);

	list.AddParameter(fa, "a");
	list.AddParameter(fb, "b");
	list.AddParameter(fc, "c");
	list.AddParameter(fd, "d");
	
	/* set the description */
	list.SetDescription("f(x) = a + b x + c (1 - exp[-x/d])");	
}

void LinearExponentialT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);

	fa = list.GetParameter("a");
	fb = list.GetParameter("b");
	fc = list.GetParameter("c");
	fd = list.GetParameter("d");

	/* check */
	if (fabs(fd) < kSmall) ExceptionT::BadInputValue("LinearExponentialT::TakeParameterList",
		"bad value of d: %d", fd);
}
