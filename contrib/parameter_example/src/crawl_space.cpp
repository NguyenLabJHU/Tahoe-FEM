/* $Id: crawl_space.cpp,v 1.1.2.2 2003-05-04 22:13:39 paklein Exp $ */
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

	list.GetParameter("sump_pump", sump_pump_);
}
