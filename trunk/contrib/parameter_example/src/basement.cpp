/* $Id: basement.cpp,v 1.3 2003-08-14 01:22:43 paklein Exp $ */
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

	length_ = list.GetParameter("length");
	width_  = list.GetParameter("width");
	height_ = list.GetParameter("height");
}
