/* $Id: FBC_ControllerT.cpp,v 1.9 2003-11-04 01:37:02 paklein Exp $ */
/* created: paklein (11/17/1997) */
#include "FBC_ControllerT.h"
#include "ArrayT.h"

#include <iostream.h>

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<FBC_ControllerT>::fByteCopy = false;
DEFINE_TEMPLATE_STATIC const bool ArrayT<FBC_ControllerT*>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
FBC_ControllerT::FBC_ControllerT(FEManagerT& fe_manager, int group):
	ParameterInterfaceT("FBC_controller"),
	fFEManager(fe_manager),
	fGroup(group),
	fIntegrator(NULL)
{

}

/* destructor */
FBC_ControllerT::~FBC_ControllerT(void) { }

/* set the controller */
void FBC_ControllerT::SetController(const eIntegratorT* controller)
{
	fIntegrator = controller;
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
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2,
	AutoArrayT<const iArray2DT*>& equivalent_nodes) const
{
#pragma unused(connects_1)
#pragma unused(connects_2)
#pragma unused(equivalent_nodes)
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
