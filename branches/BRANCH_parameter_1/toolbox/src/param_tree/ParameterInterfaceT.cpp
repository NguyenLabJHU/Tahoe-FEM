/* $Id: ParameterInterfaceT.cpp,v 1.2.2.1 2003-04-27 22:13:45 paklein Exp $ */
#include "ParameterInterfaceT.h"
#include "ParameterListT.h"

using namespace Tahoe;

/* constructor */
ParameterInterfaceT::ParameterInterfaceT(const StringT& name):
	fName(name)
{

}

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

/* return the list of subordinate names */
void ParameterInterfaceT::SubNames(ArrayT<StringT>& names, ArrayT<ParameterListT::OccurrenceT>& occur,
	ArrayT<bool>& is_inline) const
{
	names.Dimension(0);
	occur.Dimension(0);
	is_inline.Dimension(0);
}
	
/* return a pointer to the ParameterInterfaceT */
ParameterInterfaceT* ParameterInterfaceT::NewSub(const StringT& list_name)
{
#pragma unused(list_name)
	return NULL;
}
