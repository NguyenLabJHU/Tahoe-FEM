/* $Id: bedroom.cpp,v 1.3 2003-08-14 01:22:43 paklein Exp $ */
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

void bedroom::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	room::TakeParameterList(list);

	/* extract parameter values */
	floor = list.GetParameter("floor");
	
	/* construct windows */
	int num_windows = list.NumLists("window");
	windows_.Dimension(num_windows);
	const ArrayT<ParameterListT>& sub_lists = list.Lists();
	num_windows = 0;
	for (int i = 0; i < sub_lists.Length(); i++)
		if (sub_lists[i].Name() == "window")
		{
			windows_[num_windows] = new window;
			windows_[num_windows]->TakeParameterList(sub_lists[i]);
			num_windows++;
		}
}

void bedroom::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	room::DefineSubs(sub_list);

	/* the window */
	sub_list.AddSub("window", ParameterListT::Any);	
}

ParameterInterfaceT* bedroom::NewSub(const StringT& list_name) const
{
	if (list_name == "window")
		return new window;
	else /* inherited */
		return room::NewSub(list_name);
}
