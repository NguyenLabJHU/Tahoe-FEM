/* $Id: basement.cpp,v 1.1.2.1 2003-05-03 09:08:27 paklein Exp $ */
#include "basement.h"

basement::basement(const StringT& name):
	ParameterInterfaceT(name)
{

}

void basement::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	LimitT bound(0, LimitT::Lower);

	ParameterT height(height_, "height");
	height.AddLimit(bound);
	list.AddParameter(height);

	ParameterT length(length_, "length");
	length.AddLimit(bound);
	list.AddParameter(length);

	ParameterT width(width_, "width");
	width.AddLimit(bound);
	list.AddParameter(width);
}

void basement::SetParameters(const ParameterListT& list)
{
#pragma unused(list)

	list.GetParameter("length", length_);
	list.GetParameter("width", width_);
	list.GetParameter("height", height_);
}
