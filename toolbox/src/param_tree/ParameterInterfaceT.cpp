/* $Id: ParameterInterfaceT.cpp,v 1.2 2003-04-26 19:10:35 paklein Exp $ */
#include "ParameterInterfaceT.h"
#include "ParameterListT.h"

using namespace Tahoe;

/* accept completed parameter list */
void ParameterInterfaceT::SetParameters(const ParameterListT& list) 
{
	const char caller[] = "ParameterInterfaceT::SetParameters";
	if (list.Name() != Name())
		ExceptionT::GeneralFail(caller, "list name \"%s\" must be \"%s\"",
			list.Name().Pointer(), Name().Pointer());
}

/* build complete parameter list description */
void ParameterInterfaceT::DefineParameters(ParameterListT& list) const
{
	const char caller[] = "ParameterInterfaceT::DefineParameters";
	if (list.Name() != Name())
		ExceptionT::GeneralFail(caller, "list name \"%s\" must be \"%s\"",
			list.Name().Pointer(), Name().Pointer());
}

/* return the list of sub-list names */
void ParameterInterfaceT::SubListNames(ArrayT<StringT>& list, ArrayT<ParameterListT::OccurrenceT>& occur) const
{
	list.Dimension(0);
	occur.Dimension(0);
}
	
/* return a pointer to the ParameterInterfaceT */
ParameterInterfaceT* ParameterInterfaceT::SubList(const StringT& list_name)
{
#pragma unused(list_name)
	return NULL;
}
