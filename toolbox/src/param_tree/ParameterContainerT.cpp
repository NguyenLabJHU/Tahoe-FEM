/* $Id: ParameterContainerT.cpp,v 1.1 2003-11-01 21:10:37 paklein Exp $ */
#include "ParameterContainerT.h"

using namespace Tahoe;

/* constructor */
ParameterContainerT::ParameterContainerT(const StringT& name):
	ParameterListT(name),
	ParameterInterfaceT(name)
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
