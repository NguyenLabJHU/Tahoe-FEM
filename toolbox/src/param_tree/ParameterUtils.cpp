/* $Id: ParameterUtils.cpp,v 1.2 2003-08-18 03:34:13 paklein Exp $ */
#include "ParameterUtils.h"

using namespace Tahoe;

/**********************************************************************
 * IntegerListT implementation
 **********************************************************************/

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
 * DoubleListT implementation
 **********************************************************************/

/* constructors */
DoubleListT::DoubleListT(const StringT& name):
	ParameterInterfaceT(name)
{

}

DoubleListT::DoubleListT(void):
	ParameterInterfaceT("DoubleList")
{

}

/* describe the parameters needed by the interface */
void DoubleListT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	/* add optional name */
	list.AddParameter(fListName, "name", ParameterListT::ZeroOrOnce);
}

/* information about subordinate parameters */
void DoubleListT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);
	
	/* zero or more integers */
	sub_list.AddSub("Double", ParameterListT::Any);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* DoubleListT::NewSub(const StringT& list_name) const
{
	if (list_name == "Double")
		return new DoubleT;
	else
		return ParameterInterfaceT::NewSub(list_name);
}

/* accept parameter list */
void DoubleListT::TakeParameterList(const ParameterListT& list)
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
	NamedValueT<int>("Integer")
{

}

IntegerT::IntegerT(const StringT& name):
	NamedValueT<int>(name)
{

}

/**********************************************************************
 * DoubleT implementation
 **********************************************************************/

/* constructors */
DoubleT::DoubleT(void):
	NamedValueT<double>("Double")
{

}

DoubleT::DoubleT(const StringT& name):
	NamedValueT<double>(name)
{

}
