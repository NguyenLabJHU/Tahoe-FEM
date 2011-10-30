/* $Id: RodMaterialT.cpp,v 1.9.6.1 2011-10-30 06:26:11 bcyansfn Exp $ */
/* created: paklein (11/20/1996) */
#include "RodMaterialT.h"

#include <iostream>
#include "ExceptionT.h"

using namespace Tahoe;

/* constructor */
RodMaterialT::RodMaterialT(double mass):
	fMass(mass)
{
	if (fMass < 0) ExceptionT::BadInputValue("RodMaterialT::RodMaterialT");
}

/* destructor */
RodMaterialT::~RodMaterialT(void)
{
	//delete fThermal;
}

#if 0
/* thermal accessors */
int RodMaterialT::ThermalScheduleNumber(void) const
{
	return fThermal->ScheduleNum();
}

void RodMaterialT::SetThermalSchedule(const ScheduleT* LTfPtr)
{
	fThermal->SetSchedule(LTfPtr);
}

/* returns true if the material has internal forces in the unloaded
* configuration, ie thermal strains */
int RodMaterialT::HasInternalStrain(void) const { return 0; }
#endif
