/* $Id: J2_C0HardeningT.cpp,v 1.1.2.2 2004-06-11 01:38:17 paklein Exp $ */
#include "J2_C0HardeningT.h"

#include "dSymMatrixT.h"
#include "iArrayT.h"
#include <math.h>

using namespace Tahoe;

/* class constants */
const double sqrt23 = sqrt(2.0/3.0);

/* constructor */
J2_C0HardeningT::J2_C0HardeningT(void):
	ParameterInterfaceT("J2_C0_hardening"),
	fIsLinear(false),
	fK(NULL)
{

}

/* destructor */
J2_C0HardeningT::~J2_C0HardeningT(void) { delete fK; };

/* information about subordinate parameter lists */
void J2_C0HardeningT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* hardening function */
	sub_list.AddSub("hardening_function_choice", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void J2_C0HardeningT::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
	SubListT& sub_sub_list) const
{
	if (sub == "hardening_function_choice")
	{
		order = ParameterListT::Choice;
	
		/* function types */
		sub_sub_list.AddSub("linear_function");
		sub_sub_list.AddSub("cubic_spline");
		sub_sub_list.AddSub("linear_exponential");
		sub_sub_list.AddSub("power_law");
	}
	else /* inherited */
		ParameterInterfaceT::DefineInlineSub(sub, order, sub_sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* J2_C0HardeningT::NewSub(const StringT& list_name) const
{
	/* try to construct C1 function */
	C1FunctionT* function = C1FunctionT::New(list_name);
	if (function)
		return function;
	else /* inherited */
		return ParameterInterfaceT::NewSub(list_name);
}

/* accept parameter list */
void J2_C0HardeningT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "J2_C0HardeningT::TakeParameterList";

	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* construct hardening function */
	const ParameterListT* hardening = list.ResolveListChoice(*this, "hardening_function_choice");
	if (hardening) {
		fK = C1FunctionT::New(hardening->Name());
		if (!fK) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", hardening->Name().Pointer());
		fK->TakeParameterList(*hardening);

		/* set flag */
		if (hardening->Name() == "linear_function") 
			fIsLinear = true;
	}
	else
		ExceptionT::GeneralFail(caller, "could not resolve \"hardening_function_choice\"");
}

double J2_C0HardeningT::YieldCondition(const dSymMatrixT& relstress, double alpha) const 
{
	return sqrt(relstress.ScalarProduct()) - sqrt23*K(alpha);
}
