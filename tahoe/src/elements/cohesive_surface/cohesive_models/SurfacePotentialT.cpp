/* $Id: SurfacePotentialT.cpp,v 1.3 2001-10-11 00:53:41 paklein Exp $ */
/* created: paklein (06/20/1999) */

#include "SurfacePotentialT.h"

/* constructor */
SurfacePotentialT::SurfacePotentialT(int ndof):
	fTraction(ndof),
	fStiffness(ndof)
{

}

/* destructor */
SurfacePotentialT::~SurfacePotentialT(void) { }

/* initialize the state variable array */
void SurfacePotentialT::InitStateVariables(dArrayT& state)
{
	int num_state = NumStateVariables();
	if (state.Length() != num_state) {
		cout << "\n SurfacePotentialT::InitStateVariables: expecting state variable array\n"
		     <<   "     length " << num_state << ", found length " << state.Length() << endl;
		throw eSizeMismatch;
	}

	/* clear */
	if (num_state > 0) state = 0.0;
}

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

void SurfacePotentialT::ComputeOutput(const dArrayT& jump, const dArrayT& state, 
	dArrayT& output)
{
#pragma unused(jump)
#pragma unused(state)
#pragma unused(output)
}

/*************************************************************************
* Protected
*************************************************************************/

/* return true if the potential has compatible (type and sequence)
* nodal output - FALSE by default */
bool SurfacePotentialT::CompatibleOutput(const SurfacePotentialT&) const
{
	return false;
}
