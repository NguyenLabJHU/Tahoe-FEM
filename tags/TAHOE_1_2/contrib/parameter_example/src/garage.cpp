/* $Id: garage.cpp,v 1.2 2003-05-04 22:49:50 paklein Exp $ */
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

	list.GetParameter("opener", opener_);
	list.GetParameter("length", length_);
	list.GetParameter("width", width_);

	const ParameterListT* window_params = list.List("window");
	if (window_params) {
		window_ = new window;
		window_->TakeParameterList(*window_params);
	}
}

void garage::SubNames(ArrayT<StringT>& names, ArrayT<ParameterListT::OccurrenceT>& occur,
		ArrayT<bool>& is_inline) const
{
	/* temporaries */
	AutoArrayT<StringT> names_tmp;
	AutoArrayT<ParameterListT::OccurrenceT> occur_tmp;
	AutoArrayT<bool> is_inline_tmp;
	
	/* inherited */
	ParameterInterfaceT::SubNames(names_tmp, occur_tmp, is_inline_tmp);

	/* the window */
	names_tmp.Append("window");
	occur_tmp.Append(ParameterListT::ZeroOrOnce);
	is_inline_tmp.Append(false);
	
	/* copy to return values */
	names = names_tmp;
	occur = occur_tmp;
	is_inline = is_inline_tmp;
}

ParameterInterfaceT* garage::NewSub(const StringT& list_name) const
{
	if (list_name == "window")
		return new window;
	else /* inherited */
		return ParameterInterfaceT::NewSub(list_name);
}
