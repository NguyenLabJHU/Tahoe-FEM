/* $Id: ParameterContainerT.cpp,v 1.2 2004-01-21 17:05:40 paklein Exp $ */
#include "ParameterContainerT.h"

using namespace Tahoe;

/* constructor */
ParameterContainerT::ParameterContainerT(const StringT& name):
	ParameterListT(name),
	ParameterInterfaceT(name),
	fSubSource(NULL)
{

}

/* define name */
void ParameterContainerT::SetName(const StringT& name)
{
	/* inherited */
	ParameterListT::SetName(name);
	ParameterInterfaceT::SetName(name);
}

void ParameterContainerT::AddSub(const StringT& name, ParameterListT::OccurrenceT occur, 
	bool is_inline)
{
	fSubs.AddSub(name, occur, is_inline);
}

void ParameterContainerT::AddSub(const SubListDescriptionT& sub)
{
	fSubs.AddSub(sub);
}

/* set source for subs not defined by the container */
void ParameterContainerT::SetSubSource(const ParameterInterfaceT* sub_source) 
{
	if (sub_source == this)
		ExceptionT::GeneralFail("ParameterContainerT::SetSubSource", 
			"the sub source for %s cannot be self", Name().Pointer());
	fSubSource = sub_source; 
};

/* a pointer to the ParameterInterfaceT of the given subordinate*/
ParameterInterfaceT* ParameterContainerT::NewSub(const StringT& list_name) const
{
	/* inherited (get from self) */
	ParameterInterfaceT* sub = ParameterInterfaceT::NewSub(list_name);
	
	/* get from sub source */
	if (!sub && fSubSource)
		sub = fSubSource->NewSub(list_name);
		
	return sub;
}

/* return the description of the given inline subordinate parameter list. */
void ParameterContainerT::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
	SubListT& sub_sub_list) const
{
	/* inherited (get from self) */
	ParameterInterfaceT::DefineInlineSub(sub, order, sub_sub_list);

	/* if list is empty, try sub source */
	if (fSubSource && sub_sub_list.Length() == 0)
		fSubSource->DefineInlineSub(sub, order, sub_sub_list);
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* describe the parameters needed by the interface. */
void ParameterContainerT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	/* register all parameters */
	for (int i = 0; i < fParameters.Length(); i++)
		list.AddParameter(fParameters[i], fParametersOccur[i]);
}

/* information about subordinate parameter lists */
void ParameterContainerT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* register all sublists */
	for (int i = 0; i < fSubs.Length(); i++)
		sub_list.AddSub(fSubs[i]);
}
