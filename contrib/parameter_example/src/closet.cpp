/* $Id: closet.cpp,v 1.1.2.2 2003-05-03 09:08:27 paklein Exp $ */
#include "closet.h"

closet::closet(void):
	room("closet"),
	has_shelf_(false),
	has_bar_(false)
{

}

void closet::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	room::DefineParameters(list);
	
	LimitT t("true", 1);
	LimitT f("false", 0);
	
	ParameterT has_shelf(has_shelf_, "has_shelf");
	has_shelf.SetDefault(false);
	list.AddParameter(has_shelf);

	ParameterT has_bar(has_bar_, "has_bar");
	has_bar.SetDefault(true);
	list.AddParameter(has_bar);
}

void closet::SetParameters(const ParameterListT& list)
{
	/* inherited */
	room::SetParameters(list);
	
	list.GetParameter("has_shelf", has_shelf_);
	list.GetParameter("has_bar", has_bar_);
}
