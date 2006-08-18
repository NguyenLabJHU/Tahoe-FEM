/* $Id: ScaledSinh.cpp,v 1.1 2006-08-18 18:43:15 thao Exp $ */

#include "ScaledSinh.h"
#include <iostream.h>
#include <math.h>
#include "ExceptionT.h"
#include "dArrayT.h"

static const double Pi = 4.0*atan(1.0);

/* constructors */

using namespace Tahoe;

ScaledSinh::ScaledSinh(double A, double B)
{ 
	if (A > kSmall)
		fieta0 = 1.0/A;
	else ExceptionT::GeneralFail("ScaledSinh::ScaledSinh", "expecting nonzero viscosity");
	SetName("scaled_sinh");

	ftau0 = B;
}

ScaledSinh::ScaledSinh(void):fieta0(0.0), ftau0(0.0) 
{ 
	SetName("scaled_sinh");
}

/* I/O */
void ScaledSinh::Print(ostream& out) const
{

	/* parameters */
	out << " inverse viscosity. . . . . . . . . . . . . . . . . . . . . = " << fieta0 << '\n';
	out << " activation stress. . . . . . . . . . . . . . . . . . . . . . = " << ftau0 << '\n';
}

void ScaledSinh::PrintName(ostream& out) const
{
	out << "    ieta0 tau0/tau Sinh(tau/tau0) \n";
}

/*
* Returning values
*/
double ScaledSinh::Function(double tau) const
{
	double x = tau/ftau0;
	
	if (x > kSmall)
		return (fieta0 * sinh(x)/x);
	else return(fieta0);
}

double ScaledSinh::DFunction(double tau) const
{
	double x = tau/ftau0;
	
	if (x > kSmall)
		return (fieta0*(tau*cosh(x) + ftau0*sinh(x))/(x*x));
	else return(0.0);
}

double ScaledSinh::DDFunction(double tau) const
{
	double x = tau/ftau0;
	
	if (x > kSmall)
	{
		double A = (tau*tau + 2.0*ftau0*ftau0)*sinh(x);
		double B  = (-2.0*tau*ftau0)*cosh(x);
		return (fieta0*(A+B)/(tau*tau*tau*ftau0));
	}
	else return(fieta0/(3.0*ftau0*ftau0));
}

/* returning values in groups */
dArrayT& ScaledSinh::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double tau = (*pr++);
		double x = tau/ftau0;
	
		if (x > kSmall)
			*pU++ = fieta0 * sinh(x)/x;
		else *pU++ = fieta0;
	}
	return(out);
}

dArrayT& ScaledSinh::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double tau = (*pr++);
		double x = tau/ftau0;
	
		if (x > kSmall)
			*pdU++ = fieta0*(tau*cosh(x) + ftau0*sinh(x))/(x*x);
		else *pdU++ = 0.0;
	}
	return(out);
}

dArrayT& ScaledSinh::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double tau = (*pr++);
		double x = tau/ftau0;

		if (x > kSmall)
		{
			double A = (tau*tau + 2.0*ftau0*ftau0)*sinh(x);
			double B  = (-2.0*tau*ftau0)*cosh(x);
			*pddU++ = (fieta0*(A+B)/(tau*tau*tau*ftau0));
		}
		else *pddU++ = fieta0/(3.0*ftau0*ftau0);
	}
	return(out);
}

void ScaledSinh::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);

	list.AddParameter(fieta0, "inverse_viscosity");
	list.AddParameter(ftau0, "activation_stress");
	
	/* set the description */
	list.SetDescription("f(tau) = ieta_0* tau0/tau * Sinh(tau/tau0)");	
}

void ScaledSinh::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);

	/* bound */
	LimitT lower(0.0, LimitT::Lower);

	fieta0 = list.GetParameter("inverse_viscosity");
	ftau0 = list.GetParameter("activation_stress");

	/* check */
	if (fieta0 < kSmall) ExceptionT::BadInputValue("ScaledSinh::TakeParameterList",
		"expecting a positive value for the inverse viscosity: %d", fieta0);
	if (ftau0 < kSmall) ExceptionT::BadInputValue("ScaledSinh::TakeParameterList",
		"expecting a positive value for the activation stress: %d", ftau0);
}

