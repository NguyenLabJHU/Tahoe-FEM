/* $Id: garage.cpp,v 1.3 2003-08-14 01:22:43 paklein Exp $ */
#include "garage.h"
#include "window.h"

garage::garage(void):
	ParameterInterfaceT("garage"),
	opener_(false),
	length_(0.0),
	width_(0.0),
	window_(NULL)
{

}

garage::~garage(void)
{
	delete window_;
}

void garage::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	ParameterT opener(opener_, "opener");
	opener.SetDefault(true);
	list.AddParameter(opener);

	LimitT bound(0, LimitT::Lower);

	ParameterT length(length_, "length");
	length.AddLimit(bound);
	length.SetDefault(15.0);
	list.AddParameter(length);

	ParameterT width(width_, "width");
	width.AddLimit(bound);
	width.SetDefault(12.0);
	list.AddParameter(width);
}

void garage::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	opener_ = list.GetParameter("opener");
	length_ = list.GetParameter("length");
	width_  = list.GetParameter("width");

	const ParameterListT* window_params = list.List("window");
	if (window_params) {
		window_ = new window;
		window_->TakeParameterList(*window_params);
	}
}

void garage::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* the window */
	sub_list.AddSub("window", ParameterListT::ZeroOrOnce);
}

ParameterInterfaceT* garage::NewSub(const StringT& list_name) const
{
	if (list_name == "window")
		return new window;
	else /* inherited */
		return ParameterInterfaceT::NewSub(list_name);
}
