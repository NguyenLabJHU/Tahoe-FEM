/* $Id: basement.cpp,v 1.2 2003-05-04 22:49:50 paklein Exp $ */
#include "basement.h"

basement::basement(const StringT& name):
	ParameterInterfaceT(name),
	height_(0.0),
	length_(0.0),
	width_(0.0)
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

void basement::TakeParameterList(const ParameterListT& list)
{
#pragma unused(list)

	list.GetParameter("length", length_);
	list.GetParameter("width", width_);
	list.GetParameter("height", height_);
}
