/* $Id: crawl_space.cpp,v 1.1.2.1 2003-05-03 09:08:27 paklein Exp $ */
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

void crawl_space::SetParameters(const ParameterListT& list)
{
	/* inherited */
	basement::SetParameters(list);

	list.GetParameter("sump_pump", sump_pump_);
}
