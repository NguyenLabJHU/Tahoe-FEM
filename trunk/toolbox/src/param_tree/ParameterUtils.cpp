/* $Id: ParameterUtils.cpp,v 1.1 2003-08-14 01:22:03 paklein Exp $ */
#include "ParameterUtils.h"

using namespace Tahoe;

/* constructors */
IntegerListT::IntegerListT(const StringT& name):
	ParameterInterfaceT(name)
{

}

IntegerListT::IntegerListT(void):
	ParameterInterfaceT("IntegerList")
{

}

/* describe the parameters needed by the interface */
void IntegerListT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	/* add optional name */
	list.AddParameter(fListName, "name", ParameterListT::ZeroOrOnce);
}

/* information about subordinate parameters */
void IntegerListT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);
	
	/* zero or more integers */
	sub_list.AddSub("Integer", ParameterListT::Any);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* IntegerListT::NewSub(const StringT& list_name) const
{
	if (list_name == "Integer")
		return new IntegerT;
	else
		return ParameterInterfaceT::NewSub(list_name);
}

/* accept parameter list */
void IntegerListT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* get name if defined */
	const ParameterT* list_name = list.Parameter("name");
	if (list_name) 
		fListName = *list_name;
}

/**********************************************************************
 * IntegerT implementation
 **********************************************************************/

/* constructors */
IntegerT::IntegerT(void):
	ParameterInterfaceT("Integer")
{

}

IntegerT::IntegerT(const StringT& name):
	ParameterInterfaceT(name)
{

}

void IntegerT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);
	
	/* value name */
	list.AddParameter(fValueName, "name", ParameterListT::ZeroOrOnce);
	list.AddParameter(fValue, "value");
}

void IntegerT::TakeParameterList(const ParameterListT& list)
{
	/* get name */
	const ParameterT* value_name = list.Parameter("name");
	if (value_name)
		fValueName = *value_name;

	/* the value */
	list.GetParameter("value", fValue);
}
