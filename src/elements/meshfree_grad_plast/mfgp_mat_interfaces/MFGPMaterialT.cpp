/* $Id: MFGPMaterialT.cpp,v 1.2 2005-04-27 20:34:14 kyonten Exp $ */
#include "MFGPMaterialT.h"
#include "MFGPMatSupportT.h"
#include "ArrayT.h"
#include "StringT.h"

#include "dSymMatrixT.h"
#include "LocalArrayT.h"
#include "ParameterContainerT.h"

namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<MFGPMaterialT>::fByteCopy = false;
DEFINE_TEMPLATE_STATIC const bool ArrayT<MFGPMaterialT*>::fByteCopy = true;
} /* namespace Tahoe */

using namespace Tahoe;

/* constructor */
MFGPMaterialT::MFGPMaterialT(void):
	ParameterInterfaceT("mfgp_material"),
	fMFGPMatSupport(NULL),
	fDensity(0.0),
	fConstraint(kNoConstraint),
	fNumDOF(0),
	fNumSD(0),
	fNumIP(0)
{

}

/* set the material support or pass NULL to clear */
void MFGPMaterialT::SetMFGPMatSupport(const MFGPMatSupportT* support)
{
	fMFGPMatSupport = support;
	if (fMFGPMatSupport) {
		fNumDOF = fMFGPMatSupport->NumDOF();
		fNumSD = fMFGPMatSupport->NumSD();
		fNumIP = fMFGPMatSupport->NumIP();
		fDensity = 0.0;
		fConstraint = kNoConstraint;
	}
	else {
		fNumDOF = 0;
		fNumSD = 0;
		fNumIP = 0;
		fDensity = 0.0;
		fConstraint = kNoConstraint;
	}
}

/* destructor */
MFGPMaterialT::~MFGPMaterialT(void) { }

/* number of element nodes in the host element group */
int MFGPMaterialT::NumElementNodes() const {
	ElementCardT& card = ElementCard(0);
	return card.NodesU().Length();
}

/* element card data */
int MFGPMaterialT::NumElements(void) const
{
	return MFGPMatSupport().NumElements();
}

int MFGPMaterialT::CurrElementNumber(void) const
{
	return MFGPMatSupport().CurrElementNumber();
}

ElementCardT& MFGPMaterialT::ElementCard(int card) const
{
	ElementCardT* the_card = MFGPMatSupport().ElementCard(card);
	if (!the_card) ExceptionT::GeneralFail("MFGPMaterialT::ElementCard");
	return *the_card;
}

ElementCardT& MFGPMaterialT::CurrentElement(void) const
{
	ElementCardT* the_card = MFGPMatSupport().CurrentElement();
	if (!the_card) ExceptionT::GeneralFail("MFGPMaterialT::CurrentElement");
	return *the_card;
}

/* storage initialization */
bool MFGPMaterialT::NeedsPointInitialization(void) const { return false; }
void MFGPMaterialT::PointInitialize(void) { /* nothing to do */ }

/* form of tangent matrix */
GlobalT::SystemTypeT MFGPMaterialT::TangentType(void) const
{
	/* symmetric by default */
	return GlobalT::kSymmetric;
}

/* apply pre-conditions at the current time step */
void MFGPMaterialT::InitStep(void) { }

/* finalize the current time step */
void MFGPMaterialT::CloseStep(void) { }

/* update/reset internal variables */
void MFGPMaterialT::UpdateHistory(void) { }
void MFGPMaterialT::ResetHistory(void) { }

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int MFGPMaterialT::NumOutputVariables(void) const { return 0; }
void MFGPMaterialT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Free();	
}
void MFGPMaterialT::ComputeOutput(dArrayT& output)
{
#pragma unused(output)
}

/* returns true if two materials have compatible output variables */
bool MFGPMaterialT::CompatibleOutput(const MFGPMaterialT& m1, 
	const MFGPMaterialT& m2)
{
	/* number of variables */
	if (m1.NumOutputVariables() != m2.NumOutputVariables())
		return false;
	/* labels */
	else
	{
		ArrayT<StringT> labels1, labels2;
		m1.OutputLabels(labels1);
		m2.OutputLabels(labels2);
		for (int i = 0; i < labels1.Length(); i++)
			if (labels1[i] != labels2[i])
				return false;

		/* compatible if execution false through */
		return true;
	}
}

/* describe the parameters needed by the interface */
void MFGPMaterialT::DefineParameters(ParameterListT& list) const
{
	/* density */
	ParameterT density(fDensity, "density");
	density.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(density);

	/* 2D constraint option */
	ParameterT constraint(ParameterT::Enumeration, "constraint_2D");
	constraint.AddEnumeration("none", kNoConstraint);
	constraint.AddEnumeration("plane_stress", kPlaneStress);
	constraint.AddEnumeration("plane_strain", kPlaneStrain);
	constraint.SetDefault(fConstraint);
	list.AddParameter(constraint);
}

/* accept parameter list */
void MFGPMaterialT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "MFGPMaterialT::TakeParameterList";

	/* density */
	fDensity = list.GetParameter("density");

	/* 2D constraint - default to plane strain for 2D materials */
	list.GetParameter("constraint_2D", enum2int<ConstraintT>(fConstraint));
	if (NumSD() == 3)
		fConstraint = kNoConstraint;
	else if (NumSD() == 2 && fConstraint == kNoConstraint)
		fConstraint = kPlaneStrain;
}
	
