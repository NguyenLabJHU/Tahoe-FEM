/* $Id: RodMaterialT.cpp,v 1.8 2004-07-15 08:30:22 paklein Exp $ */
/* created: paklein (11/20/1996) */
#include "RodMaterialT.h"

#include <iostream.h>

#include "ExceptionT.h"
#include "ThermalDilatationT.h"
#include "ifstreamT.h"

using namespace Tahoe;

#pragma message("XML clean up")

/* constructor */
RodMaterialT::RodMaterialT(ifstreamT& in)
{
	fMass = -1.0;
	in >> fMass;
	if (fMass < 0) throw ExceptionT::kBadInputValue;
	
//	fThermal = new ThermalDilatationT(in);
	fThermal = new ThermalDilatationT;
	if (!fThermal) throw ExceptionT::kOutOfMemory;
}

/* destructor */
RodMaterialT::~RodMaterialT(void)
{
	delete fThermal;
}

/* print parameters */
void RodMaterialT::PrintParameters(ostream& out) const
{
	/* thermal expansion parameters */
//	fThermal->Print(out);
}

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
