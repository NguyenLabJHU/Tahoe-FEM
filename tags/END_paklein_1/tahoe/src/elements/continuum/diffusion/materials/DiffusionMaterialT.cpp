/* $Id: DiffusionMaterialT.cpp,v 1.3.4.1 2002-10-17 04:38:19 paklein Exp $ */
/* created: paklein (10/02/1999)                                          */

#include "DiffusionMaterialT.h"

#include <iostream.h>

#include "StringT.h"
#include "fstreamT.h"
#include "dArrayT.h"
#include "dSymMatrixT.h"

#include "LocalArrayT.h"
#include "DiffusionT.h"

/* constructor */

using namespace Tahoe;

DiffusionMaterialT::DiffusionMaterialT(ifstreamT& in, const DiffusionT& element):
	ContinuumMaterialT(element),
	fLocDisp(element.Displacements()),
	fConductivity(NumSD()),
	fT_x(1, NumSD()),
	fq_i(NumSD())
{
	in >> fDensity;		 if (fDensity <= 0.0) throw ExceptionT::kBadInputValue;
	in >> fSpecificHeat; if (fDensity <= 0.0) throw ExceptionT::kBadInputValue;
	in >> fConductivity;

	fCapacity = fDensity*fSpecificHeat;
}

/* I/O functions */
void DiffusionMaterialT::Print(ostream& out) const
{
	/* inherited */
	ContinuumMaterialT::Print(out);

	out << " Density . . . . . . . . . . . . . . . . . . . . = " << fDensity      << '\n';
	out << " Specific Heat . . . . . . . . . . . . . . . . . = " << fSpecificHeat << '\n';
	out << " Conductivity:\n" << fConductivity << endl;
}

/* heat flux */
const dArrayT& DiffusionMaterialT::q_i(void)
{
	/* should be 1 row */
	ContinuumElement().IP_ComputeGradient(fLocDisp, fT_x);
	fConductivity.Multx(fT_x, fq_i);
	return fq_i;
}

/*************************************************************************
* Protected
*************************************************************************/

void DiffusionMaterialT::PrintName(ostream& out) const
{
	/* inherited */
	ContinuumMaterialT::PrintName(out);
	
	out << "    Linear diffusion material\n";
}
