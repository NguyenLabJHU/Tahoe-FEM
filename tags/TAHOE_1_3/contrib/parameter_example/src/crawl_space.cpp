/* $Id: crawl_space.cpp,v 1.3 2003-08-14 01:22:43 paklein Exp $ */
#include "crawl_space.h"

crawl_space::crawl_space(void):
	basement("crawl_space"),
	sump_pump_(true)
{

}

void crawl_space::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	basement::DefineParameters(list);

	ParameterT sump_pump(sump_pump_, "sump_pump");
	sump_pump.SetDefault(true);
	list.AddParameter(sump_pump);
}

void crawl_space::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	basement::TakeParameterList(list);

	sump_pump_ = list.GetParameter("sump_pump");
}
