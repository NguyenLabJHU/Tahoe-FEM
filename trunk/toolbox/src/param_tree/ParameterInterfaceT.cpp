/* $Id: ParameterInterfaceT.cpp,v 1.1 2003-04-22 18:32:16 paklein Exp $ */
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
void ParameterInterfaceT::DefineParameters(ParameterListT& list) 
{
	const char caller[] = "ParameterInterfaceT::DefineParameters";
	if (list.Name() != Name())
		ExceptionT::GeneralFail(caller, "list name \"%s\" must be \"%s\"",
			list.Name().Pointer(), Name().Pointer());
}
