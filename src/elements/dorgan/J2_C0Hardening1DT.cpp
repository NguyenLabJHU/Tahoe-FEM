/* $Id: J2_C0Hardening1DT.cpp,v 1.2 2004-08-04 22:02:14 rdorgan Exp $ */
#include "J2_C0Hardening1DT.h"

#include "dSymMatrixT.h"
#include "iArrayT.h"
#include <math.h>

using namespace Tahoe;

/* class constants */
const double sqrt23 = sqrt(2.0/3.0);

/* constructor */
J2_C0Hardening1DT::J2_C0Hardening1DT(void):
	ParameterInterfaceT("J2_C0_hardening_1D"),
	fIsLinear(false),
	fK(NULL)
{

}

/* destructor */
J2_C0Hardening1DT::~J2_C0Hardening1DT(void) { delete fK; };

/* information about subordinate parameter lists */
void J2_C0Hardening1DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* hardening function */
	sub_list.AddSub("hardening_function_choice", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void J2_C0Hardening1DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "hardening_function_choice")
	{
		order = ParameterListT::Choice;
	
		/* function types */
		sub_lists.AddSub("linear_function");
		sub_lists.AddSub("cubic_spline");
		sub_lists.AddSub("linear_exponential");
		sub_lists.AddSub("power_law");
		sub_lists.AddSub("piecewise_linear");
	}
	else /* inherited */
		ParameterInterfaceT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* J2_C0Hardening1DT::NewSub(const StringT& name) const
{
	/* try to construct C1 function */
	C1FunctionT* function = C1FunctionT::New(name);
	if (function)
		return function;
	else /* inherited */
		return ParameterInterfaceT::NewSub(name);
}

/* accept parameter list */
void J2_C0Hardening1DT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "J2_C0Hardening1DT::TakeParameterList";

	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* construct hardening function */
	const ParameterListT& hardening = list.GetListChoice(*this, "hardening_function_choice");
	fK = C1FunctionT::New(hardening.Name());
	if (!fK) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", hardening.Name().Pointer());
	fK->TakeParameterList(hardening);

	/* set flag */
	if (hardening.Name() == "linear_function") 
		fIsLinear = true;
}

double J2_C0Hardening1DT::YieldCondition(const dSymMatrixT& relstress, double alpha) const 
{
	return sqrt(relstress.ScalarProduct()) - K(alpha);
}
