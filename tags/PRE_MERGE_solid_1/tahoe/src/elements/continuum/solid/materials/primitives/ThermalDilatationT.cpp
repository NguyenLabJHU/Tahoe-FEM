/* $Id: ThermalDilatationT.cpp,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (08/25/1996)                                          */

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

/* set LTf pointer */
int ThermalDilatationT::LTfNumber(void) const { return LTfnum; }
void ThermalDilatationT::SetLTfPtr(const LoadTime* LTf) { LTfPtr = LTf; }

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
