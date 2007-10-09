/* $Id: PowerTrig.cpp,v 1.2 2007-10-09 23:24:45 rjones Exp $ */

#include "PowerTrig.h"
//#include <iostream.h>
#include <math.h>
#include "ExceptionT.h"
#include "dArrayT.h"

static const double Pi = 4.0*atan(1.0);

/* constructors */

using namespace Tahoe;

PowerTrig::PowerTrig(double A, double B, double C, double n, double phi): 
	fA(A), 
	fB(B),
	fC(C),
	fn(n),
	fphi(phi) 
{ 
	SetName("power_trig");
}

PowerTrig::PowerTrig(void): 
	fA(0.0), 
	fB(0.0),
	fC(0.0),
	fn(0.0),
	fphi(0.0) 
{ 
	SetName("power_trig");
}

/*
* Returning values
*/
double PowerTrig::Function(double r) const
{
	double cost = cos(r+fphi);
	double sint = sin(r+fphi);
	double cospt = pow(cost,fn);
	double sinpt = pow(sint,fn);
	return (fA*cospt+fB*sinpt+fC);
}

double PowerTrig::DFunction(double r) const
{
	double cost = cos(r+fphi);
	double sint = sin(r+fphi);
	double cospt = pow(cost,(fn-1.0));
	double sinpt = pow(sint,(fn-1.0));
	return(-fA*fn*cospt*sint+fB*fn*cost*sinpt);
}

double PowerTrig::DDFunction(double r) const
{
	double cost = cos(r+fphi);
	double sint = sin(r+fphi);
	double cospt = pow(cost,(fn-2.0));
	double sinpt = pow(sint,(fn-2.0));
	
	double val = fA*(fn-1.0)*fn*cospt*sint*sint+fB*(fn-1.0)*fn*sinpt*cost*cost;
	
	return(val-fA*fn*cospt*cost*cost-fB*fn*sinpt*sint*sint);
}

void PowerTrig::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);

	list.AddParameter(fA, "a");
	list.AddParameter(fB, "b");
	list.AddParameter(fC, "c");
	list.AddParameter(fn, "n");
	list.AddParameter(fphi, "phi");
	
	/* set the description */
	list.SetDescription("f(theta) = a*cos(theta+phi)^n + b*sin(theta+phi)^n + c\n");	
}

void PowerTrig::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);

	fA = list.GetParameter("a");
	fB = list.GetParameter("b");
	fC = list.GetParameter("c");
	fn = list.GetParameter("n");
	fphi = list.GetParameter("phi");
	/* check */
	if (fn < kSmall) ExceptionT::BadInputValue("PowerTrig::TakeParameterList",
		"bad value of d: %d", fn);
}
