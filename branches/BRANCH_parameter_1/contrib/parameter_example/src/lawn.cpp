/* $Id: lawn.cpp,v 1.1.2.2 2003-05-03 09:08:27 paklein Exp $ */
#include "lawn.h"

lawn::lawn(void):
	ParameterInterfaceT("lawn")
{

}

void lawn::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);
}

void lawn::SetParameters(const ParameterListT& list)
{
#pragma unused(list)
}
