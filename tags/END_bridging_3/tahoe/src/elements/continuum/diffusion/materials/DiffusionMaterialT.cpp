/* $Id: DiffusionMaterialT.cpp,v 1.6 2003-06-09 06:54:30 paklein Exp $ */
/* created: paklein (10/02/1999) */
#include "DiffusionMaterialT.h"
#include "DiffusionMatSupportT.h"

#include "StringT.h"
#include "fstreamT.h"
#include "dArrayT.h"
#include "dSymMatrixT.h"

using namespace Tahoe;

/* constructor */
DiffusionMaterialT::DiffusionMaterialT(ifstreamT& in, const DiffusionMatSupportT& support):
	ContinuumMaterialT(support),
	fDiffusionMatSupport(support),
	fConductivity(NumSD()),
	fq_i(NumSD())
{
	in >> fDensity;		 if (fDensity <= 0.0) throw ExceptionT::kBadInputValue;
	in >> fSpecificHeat; if (fDensity <= 0.0) throw ExceptionT::kBadInputValue;
	in >> fConductivity;
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
	fConductivity.Multx(fDiffusionMatSupport.Gradient(), fq_i, -1.0);
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