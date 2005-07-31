/* $Id: DiffusionMaterialT.cpp,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (10/02/1999)                                          */

#include "DiffusionMaterialT.h"

#include <iostream.h>

#include "StringT.h"
#include "fstreamT.h"
#include "dArrayT.h"
#include "dSymMatrixT.h"
#include "ShapeFunctionT.h"
#include "LocalArrayT.h"
#include "DiffusionT.h"

/* constructor */
DiffusionMaterialT::DiffusionMaterialT(ifstreamT& in, const DiffusionT& element):
	ContinuumMaterialT(element),
	fShapes(element.ShapeFunction()),
	fLocDisp(element.Displacements()),
	fConductivity(fShapes.NumSD()),
	fT_x(1,fShapes.NumSD()),
	fq_i(fShapes.NumSD())
{
	in >> fDensity;		 if (fDensity <= 0.0) throw eBadInputValue;
	in >> fSpecificHeat; if (fDensity <= 0.0) throw eBadInputValue;
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
	fShapes.GradU(fLocDisp, fT_x);
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
