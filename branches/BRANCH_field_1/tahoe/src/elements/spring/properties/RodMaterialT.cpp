/* $Id: RodMaterialT.cpp,v 1.3 2001-12-17 00:12:01 paklein Exp $ */
/* created: paklein (11/20/1996) */

#include "RodMaterialT.h"
#include <iostream.h>
#include "ExceptionCodes.h"
#include "ThermalDilatationT.h"

/* constructor */
RodMaterialT::RodMaterialT(ifstreamT& in)
{
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
	fThermal->Print(out);	
}

/* print parameters */
void RodMaterialT::PrintParameters(ostream& out) const
{
	/* thermal expansion parameters */
	fThermal->Print(out);
}

/* thermal accessors */
int RodMaterialT::ThermalLTfNumber(void) const
{
	return fThermal->LTfNumber();
}

void RodMaterialT::SetThermalLTfPtr(const LoadTime* LTfPtr)
{
	fThermal->SetLTfPtr(LTfPtr);
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
