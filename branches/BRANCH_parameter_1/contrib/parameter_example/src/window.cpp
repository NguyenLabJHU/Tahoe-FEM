/* $Id: window.cpp,v 1.1.2.2 2003-05-03 17:45:23 paklein Exp $ */
#include "window.h"

window::window(void):
  ParameterInterfaceT("window"),
	width_(0.0),
	height_(0.0)
{

}

void window::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	LimitT bound(0, LimitT::Lower);

	ParameterT width(width_, "width");
	width.AddLimit(bound);
	width.SetDefault(2.5);
	list.AddParameter(width);

	ParameterT height(height_, "height");
	height.AddLimit(bound);
	height.SetDefault(4.0);
	list.AddParameter(height);
}

void window::SetParameters(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::SetParameters(list);

	list.GetParameter("height", height_);
	list.GetParameter("width", width_);
}
