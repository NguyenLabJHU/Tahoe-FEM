/* $Id: RodMaterialT.cpp,v 1.5 2002-06-29 16:17:22 paklein Exp $ */
/* created: paklein (11/20/1996) */

#include "RodMaterialT.h"
#include <iostream.h>
#include "ExceptionCodes.h"
#include "ThermalDilatationT.h"
#include "ifstreamT.h"

/* constructor */
RodMaterialT::RodMaterialT(ifstreamT& in)
{
	fMass = -1.0;
	in >> fMass;
	if (fMass < 0) throw eBadInputValue;
	
	fThermal = new ThermalDilatationT(in);
	if (!fThermal) throw(eOutOfMemory);
}

/* destructor */
RodMaterialT::~RodMaterialT(void)
{
	delete fThermal;
}

/* print parameters */
void RodMaterialT::Print(ostream& out) const
{
	out << " Mass . . . . . . . . . . . . . . . . . . . . . .= " << fMass << '\n';	
	fThermal->Print(out);	
}

/* print parameters */
void RodMaterialT::PrintParameters(ostream& out) const
{
	/* thermal expansion parameters */
	fThermal->Print(out);
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

/***********************************************************************
* Protected
***********************************************************************/

void RodMaterialT::PrintName(ostream& out) const
{
	out << " Material name:\n";
}
