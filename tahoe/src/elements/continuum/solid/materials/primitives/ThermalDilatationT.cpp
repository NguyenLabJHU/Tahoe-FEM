/* $Id: ThermalDilatationT.cpp,v 1.4.56.1 2004-04-08 07:33:18 paklein Exp $ */
/* created: paklein (08/25/1996) */
#include "ThermalDilatationT.h"
#include "ScheduleT.h"

using namespace Tahoe;

/* constructor */
ThermalDilatationT::ThermalDilatationT(void):
	LTfPtr(NULL),
	fPercentElongation(0.0),
	LTfnum(-1)
{

}

/* returns the current elongation factor */
double ThermalDilatationT::PercentElongation(void) const
{
	return (LTfPtr != NULL) ?
		fPercentElongation*0.01*(LTfPtr->Value()) :
		0.0;
}
