/* $Id: SurfacePotentialT.cpp,v 1.18 2004-06-17 07:13:28 paklein Exp $ */
/* created: paklein (06/20/1999) */
#include "SurfacePotentialT.h"

using namespace Tahoe;

/* constructor */
SurfacePotentialT::SurfacePotentialT(int ndof):
	fTraction(ndof),
	fStiffness(ndof)
{

}

/* destructor */
SurfacePotentialT::~SurfacePotentialT(void) { }

/* initialize the state variable array */
void SurfacePotentialT::InitStateVariables(ArrayT<double>& state)
{
	int num_state = NumStateVariables();
	if (state.Length() != num_state) {
#ifndef _FRACTURE_INTERFACE_LIBRARY_	
		cout << "\n SurfacePotentialT::InitStateVariables: expecting state variable array\n"
		     <<   "     length " << num_state << ", found length " << state.Length() << endl;
#endif
		throw ExceptionT::kSizeMismatch;	
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

double SurfacePotentialT::IncrementalHeat(const dArrayT& jump, const ArrayT<double>& state)
{
#pragma unused(jump)
#pragma unused(state)
	return 0.0;
}

void SurfacePotentialT::OutputLabels(ArrayT<StringT>& labels) const
{
#pragma unused(labels)
}

void SurfacePotentialT::ComputeOutput(const dArrayT& jump, const ArrayT<double>& state, 
	dArrayT& output)
{
#pragma unused(jump)
#pragma unused(state)
#pragma unused(output)
}

/*bool SurfacePotentialT::NeedsNodalInfo(void)
{
        return false;
}

int SurfacePotentialT::NodalQuantityNeeded(void) 
{
	return 0;
}


double SurfacePotentialT::ComputeNodalValue(const dArrayT& nodalRow)
{
#pragma unused(nodalRow)
	return 0.0;
}

void SurfacePotentialT::UpdateStateVariables(const dArrayT& IPdata, ArrayT<double>& state)
{
#pragma unused(IPdata)
#pragma unused(state)
}


int SurfacePotentialT::ElementGroupNeeded(void)
{
	return -1;
}
*/

/*************************************************************************
 * Protected
 *************************************************************************/

/* return true if the potential has compatible (type and sequence)
* nodal output - FALSE by default */
bool SurfacePotentialT::CompatibleOutput(const SurfacePotentialT&) const
{
	return false;
}
