/* $Id: room.cpp,v 1.1.2.2 2003-05-03 09:08:27 paklein Exp $ */
#include "room.h"

room::room(const StringT& name):
	ParameterInterfaceT(name),
	length(0),
	width(0)
{

}

void room::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	LimitT bound(0, LimitT::Lower);

	ParameterT the_length(length, "length");
	the_length.AddLimit(bound);
	list.AddParameter(the_length);

	ParameterT the_width(width, "width");
	the_width.AddLimit(bound);
	list.AddParameter(the_width);
}

void room::SetParameters(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::SetParameters(list);
	
	list.GetParameter("length", length);
	list.GetParameter("width", width);
}
