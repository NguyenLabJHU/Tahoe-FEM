/* $Id: lawn.cpp,v 1.2 2003-05-04 22:49:50 paklein Exp $ */
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

void lawn::TakeParameterList(const ParameterListT& list)
{
#pragma unused(list)
}
