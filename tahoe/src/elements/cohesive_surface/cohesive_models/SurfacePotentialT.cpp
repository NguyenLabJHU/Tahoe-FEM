/* $Id: SurfacePotentialT.cpp,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (06/20/1999)                                          */
/* base class for surface potential with jump vector arguments            */

#include "SurfacePotentialT.h"

/* constructor */
SurfacePotentialT::SurfacePotentialT(int ndof):
	fTraction(ndof),
	fStiffness(ndof)
{

}

/* destructor */
SurfacePotentialT::~SurfacePotentialT(void) { }

/* returns true if two materials have compatible nodal outputs */
bool SurfacePotentialT::CompatibleOutput(const SurfacePotentialT& pot1,
	const SurfacePotentialT& pot2)
{
	return pot2.CompatibleOutput(pot1) && pot1.CompatibleOutput(pot2);
}
	
/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int SurfacePotentialT::NumOutputVariables(void) const
{
	return 0;
}

void SurfacePotentialT::OutputLabels(ArrayT<StringT>& labels) const
{
#pragma unused(labels)
}

void SurfacePotentialT::ComputeOutput(const dArrayT& jump_u, dArrayT& output)
{
#pragma unused(jump_u)
#pragma unused(output)
}

/*************************************************************************
* Protected
*************************************************************************/

/* return true if the potential has compatible (type and sequence)
* nodal output - FALSE by default */
bool SurfacePotentialT::CompatibleOutput(const SurfacePotentialT& potential) const
{
#pragma unused(potential)
	return false;
}
