/* $Id: closet.cpp,v 1.3 2003-08-14 01:22:43 paklein Exp $ */
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

void closet::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	room::TakeParameterList(list);
	
	has_shelf_ = list.GetParameter("has_shelf");
	has_bar_ = list.GetParameter("has_bar");
}
