/* $Id: FBC_ControllerT.cpp,v 1.14 2006-06-19 15:25:34 r-jones Exp $ */
/* created: paklein (11/17/1997) */
#include "FBC_ControllerT.h"
#include "ArrayT.h"
#include "FieldT.h"
#include "IntegratorT.h"

#include <iostream.h>

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<FBC_ControllerT>::fByteCopy = false;
DEFINE_TEMPLATE_STATIC const bool ArrayT<FBC_ControllerT*>::fByteCopy = true;
} /* namespace Tahoe */

/* converts strings to FBC_ControllerT::CodeT */
FBC_ControllerT::CodeT FBC_ControllerT::Code(const char* name)
{
	if (strcmp("sphere_penalty", name) == 0)
		return kPenaltySphere;
	else if (strcmp("wall_penalty", name) == 0)
		return kPenaltyWall;
	else if (strcmp("sphere_augmented_Lagrangian", name) == 0)
		return kAugLagSphere;
	else if (strcmp("wall_augmented_Lagrangian", name) == 0)
		return kAugLagWall;
	else if (strcmp("sphere_penalty_meshfree", name) == 0)
		return kMFPenaltySphere;
	else if (strcmp("cylinder_penalty", name) == 0)
		return kPenaltyCylinder;
	else if (strcmp("augmented_Lagrangian_KBC_meshfree", name) == 0)
		return kMFAugLagMult;
	else if (strcmp("cylinder_augmented_Lagrangian", name) == 0)
		return kAugLagCylinder;
	else if (strcmp("field_augmented_Lagrangian_KBC_meshfree", name) == 0)
		return kFieldMFAugLagMult;
	else if (strcmp("pressure_bc", name) == 0)
		return kPressureBC;
	else
		return kNone;
}

FBC_ControllerT::FBC_ControllerT(void):
	ParameterInterfaceT("FBC_controller"),
	fFieldSupport(NULL),
	fField(NULL),
	fGroup(-1),
	fIntegrator(NULL)
{

}

/* destructor */
FBC_ControllerT::~FBC_ControllerT(void) { }

/* set the associated field */
void FBC_ControllerT::SetField(const FieldT& field)
{
	/* the field */
	fField = &field;

	/* the solver group */
	fGroup = fField->Group();

	/* the support */
	fFieldSupport = &(fField->FieldSupport());

	/* the integrator */
	const IntegratorT& integrator = fField->Integrator();
	fIntegrator = &(integrator.eIntegrator());
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

/* returns true if the internal force has been changed since the last time step */
GlobalT::RelaxCodeT FBC_ControllerT::RelaxSystem(void) { return GlobalT::kNoRelax; }

/* accept parameter list */
void FBC_ControllerT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);	

	/* field support should already be set */
	if (!fFieldSupport)
		ExceptionT::GeneralFail("FBC_ControllerT::TakeParameterList", "call SetField first");
}
