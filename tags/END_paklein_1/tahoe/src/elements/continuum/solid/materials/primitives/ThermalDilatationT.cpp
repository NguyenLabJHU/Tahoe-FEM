/* $Id: ThermalDilatationT.cpp,v 1.4 2002-07-02 19:56:22 cjkimme Exp $ */
/* created: paklein (08/25/1996) */

#include "ThermalDilatationT.h"

#include <iostream.h>

#include "fstreamT.h"
#include "ScheduleT.h"

/* constructor */

using namespace Tahoe;

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
		fPercentElongation*0.01*(LTfPtr->Value()) :
		0.0;
}
