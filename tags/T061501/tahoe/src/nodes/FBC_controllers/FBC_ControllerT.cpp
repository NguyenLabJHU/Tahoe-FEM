/* $Id: FBC_ControllerT.cpp,v 1.1.1.1 2001-01-29 08:20:40 paklein Exp $ */
/* created: paklein (11/17/1997)                                          */
/* Base class for all force BC controllers                                */

#include "FBC_ControllerT.h"

#include <iostream.h>

/* constructor */
FBC_ControllerT::FBC_ControllerT(FEManagerT& fe_manager):
	fFEManager(fe_manager),
	fController(NULL)
{

}

/* destructor */
FBC_ControllerT::~FBC_ControllerT(void) { }

/* set the controller */
void FBC_ControllerT::SetController(eControllerT* controller)
{
	fController = controller;
}

/* append element equations numbers to the list */
void FBC_ControllerT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_1)
#pragma unused(eq_2)
// By default, the FBC controllers do not generate any additional
// degrees of freedom and therefore do not need to send any equation
// sets to the solver.
}

void FBC_ControllerT::Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
#pragma unused(connects_1)
#pragma unused(connects_2)
// By default, the FBC controllers do not generate any additional
// degrees of freedom and therefore do not need to send any DOF tag
// sets to the solver.
}

void FBC_ControllerT::ReadRestart(istream& in)
{
#pragma unused(in)
}

void FBC_ControllerT::WriteRestart(ostream& out) const
{
#pragma unused(out)
}
