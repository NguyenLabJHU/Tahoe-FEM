/* $Id: bedroom.cpp,v 1.1.2.3 2003-05-03 17:45:23 paklein Exp $ */
#include "bedroom.h"
#include "window.h"

bedroom::bedroom(void):
	room("bedroom"),
	floor(0)
{

}

bedroom::~bedroom(void)
{
	for (int i = 0; i < windows_.Length(); i++)
		delete windows_[i];
}

void bedroom::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	room::DefineParameters(list);

	list.AddParameter(floor, "floor");
}

void bedroom::SetParameters(const ParameterListT& list)
{
	/* inherited */
	room::SetParameters(list);

	/* extract parameter values */
	list.GetParameter("floor", floor);
	
	/* construct windows */
	int num_windows = list.NumLists("window");
	windows_.Dimension(num_windows);
	const ArrayT<ParameterListT>& sub_lists = list.Lists();
	num_windows = 0;
	for (int i = 0; i < sub_lists.Length(); i++)
		if (sub_lists[i].Name() == "window")
		{
			windows_[num_windows] = new window;
			windows_[num_windows]->SetParameters(sub_lists[i]);
			num_windows++;
		}
}

void bedroom::SubNames(ArrayT<StringT>& names, ArrayT<ParameterListT::OccurrenceT>& occur,
		ArrayT<bool>& is_inline) const
{
	/* temporaries */
	AutoArrayT<StringT> names_tmp;
	AutoArrayT<ParameterListT::OccurrenceT> occur_tmp;
	AutoArrayT<bool> is_inline_tmp;
	
	/* inherited */
	room::SubNames(names_tmp, occur_tmp, is_inline_tmp);

	/* the window */
	names_tmp.Append("window");
	occur_tmp.Append(ParameterListT::Any);
	is_inline_tmp.Append(false);
	
	/* copy to return values */
	names = names_tmp;
	occur = occur_tmp;
	is_inline = is_inline_tmp;
}

ParameterInterfaceT* bedroom::NewSub(const StringT& list_name) const
{
	if (list_name == "window")
		return new window;
	else /* inherited */
		return room::NewSub(list_name);
}
