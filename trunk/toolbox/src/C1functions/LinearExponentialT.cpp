/* $Id: LinearExponentialT.cpp,v 1.3 2002-10-20 22:38:48 paklein Exp $ */
/* created: paklein (10/30/1997)                                          */

#include "LinearExponentialT.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

/** constructor */

using namespace Tahoe;

LinearExponentialT::LinearExponentialT(double a, double b, double c, double d):
	fa(a),
	fb(b),
	fc(c),
	fd(d)
{
	if (fabs(fd) < kSmall) throw ExceptionT::kBadInputValue;
}

/** print parameters */
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

/** print function name */
void LinearExponentialT::PrintName(ostream& out) const
{
	out << "    Linear-exponential\n";
}

/** evaluate function */
double LinearExponentialT::Function(double x) const
{
	return fa + fb*x + fc*(1.0 - exp(-x/fd));
}

/** evaluate first derivative function */
double LinearExponentialT::DFunction(double x) const
{
	return fb + fc*exp(-x/fd)/fd;
}

/** evaluate second derivative function */
double LinearExponentialT::DDFunction(double x) const
{
	return -fc*exp(-x/fd)/fd/fd;
}

/* Returning values in groups */

/** multiple function evaluations */
dArrayT& LinearExponentialT::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	int length = in.Length();
	double* y = out.Pointer();
	double* x = in.Pointer();
	for (int i = 0; i < length; i++)
	{
		*y++ = fa + fb*(*x) + fc*(1.0 - exp(-(*x)/fd));
		x++;
	}
	return out;
}

/** multiple first derivative evaluations */
dArrayT& LinearExponentialT::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	int length = in.Length();
	double* y = out.Pointer();
	double* x = in.Pointer();
	for (int i = 0; i < length; i++)
		*y++ = fb + fc*exp(-(*x++)/fd)/fd;
	return out;
}

/** multiple second derivative evaluations */
dArrayT& LinearExponentialT::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	int length = in.Length();
	double* y = out.Pointer();
	double* x = in.Pointer();
	for (int i = 0; i < length; i++)
		*y++ = -fc*exp(-(*x++)/fd)/fd/fd;
	return out;
}
