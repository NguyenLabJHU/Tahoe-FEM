/* $Id: storm_shelter.cpp,v 1.3 2003-08-14 01:22:43 paklein Exp $ */
#include "storm_shelter.h"

storm_shelter::storm_shelter(void):
	basement("storm_shelter"),
	auxiliary_power_(false),
	first_aid_kit_(true),
	stored_water_(0.0)
{

}

void storm_shelter::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	basement::DefineParameters(list);

	ParameterT auxiliary_power(auxiliary_power_, "auxiliary_power");
	auxiliary_power.SetDefault(true);
	list.AddParameter(auxiliary_power);

	ParameterT first_aid_kit(first_aid_kit_, "first_aid_kit");
	first_aid_kit.SetDefault(true);
	list.AddParameter(first_aid_kit);

	ParameterT stored_water(stored_water_, "stored_water");
	stored_water.SetDefault(0.0);
	list.AddParameter(stored_water);
}

void storm_shelter::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	basement::TakeParameterList(list);

	auxiliary_power_ = list.GetParameter("auxiliary_power");
	first_aid_kit_ = list.GetParameter("first_aid_kit");
	stored_water_ = list.GetParameter("stored_water");
}
