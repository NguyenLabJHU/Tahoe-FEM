/* $Id: ScheduleT.cpp,v 1.5.22.1 2004-04-08 07:33:55 paklein Exp $ */
/* created: paklein (05/24/1996) */
#include "ScheduleT.h"

/* C1 functions */
#include "PiecewiseLinearT.h"

using namespace Tahoe;

/* constructor */
ScheduleT::ScheduleT(void):
	ParameterInterfaceT("schedule_function"),
	fFunction(NULL)
{

}

ScheduleT::ScheduleT(double value):
	ParameterInterfaceT("schedule_function"),
	fFunction(NULL)
{
	/* construct constant function */
	dArray2DT points(1,2);
	points(0,0) = 0.0;
	points(0,1) = value;
	fFunction = new PiecewiseLinearT(points);
}

/* destructor */
ScheduleT::~ScheduleT(void) { delete fFunction; }

/* set the load factor based on the time given */
void ScheduleT::SetTime(double time)
{
	fCurrentTime = time;
	fCurrentValue = fFunction->Function(fCurrentTime);
}

double ScheduleT::Value(double time) const { return fFunction->Function(time); }

/* information about subordinate parameter lists */
void ScheduleT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);
	
	/* C1FunctionT choice */
	sub_list.AddSub("function_choice", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void ScheduleT::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
	SubListT& sub_sub_list) const
{
	if (sub == "function_choice")
	{
		order = ParameterListT::Choice;
	
		/* function types */
		sub_sub_list.AddSub("piecewise_linear");
		sub_sub_list.AddSub("cubic_spline");
	}
	else /* inherited */
		return ParameterInterfaceT::DefineInlineSub(sub, order, sub_sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ScheduleT::NewSub(const StringT& list_name) const
{
	/* try to construct C1 function */
	C1FunctionT* function = C1FunctionT::New(list_name);
	if (function)
		return function;
	else /* inherited */
		return ParameterInterfaceT::NewSub(list_name);
}

/* accept parameter list */
void ScheduleT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* clear any previous */
	delete fFunction;
	fFunction = NULL;

	/* try to construct C1 function */
	const ArrayT<ParameterListT>& subs = list.Lists();
	for (int i = 0; i < subs.Length() && !fFunction; i++)
	{
		fFunction = C1FunctionT::New(subs[i].Name());
		if (fFunction)
			fFunction->TakeParameterList(subs[i]);
	}

	/* failed */
	if (!fFunction) ExceptionT::GeneralFail("ScheduleT::TakeParameterList");
}
