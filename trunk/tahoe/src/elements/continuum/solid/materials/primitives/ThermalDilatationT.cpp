/* $Id: ThermalDilatationT.cpp,v 1.2 2001-09-15 01:17:36 paklein Exp $ */
/* created: paklein (08/25/1996) */

#include "ThermalDilatationT.h"

#include <iostream.h>

#include "fstreamT.h"
#include "LoadTime.h"

/* constructor */
ThermalDilatationT::ThermalDilatationT(ifstreamT& in):
	LTfPtr(NULL)
{
	in >> LTfnum; LTfnum--;
	in >> fPercentElongation;
	
	/* overwrite */
	if (LTfnum == -1) fPercentElongation = 0.0;
}

/* I/O functions */
void ThermalDilatationT::Print(ostream& out) const
{
	out << " Dilatation LTf. . . . . . . . . . . . . . . . . = " << LTfnum + 1         << '\n';
	out << " Percent elongation. . . . . . . . . . . . . . . = " << fPercentElongation << '\n';
}

/* returns the current elongation factor */
double ThermalDilatationT::PercentElongation(void) const
{
	return (LTfPtr != NULL) ?
		fPercentElongation*0.01*(LTfPtr->LoadFactor()) :
		0.0;
}
