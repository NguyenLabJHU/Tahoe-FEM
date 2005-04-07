/* $Id: ContinuumMaterialT.cpp,v 1.12 2004-10-21 18:50:04 paklein Exp $ */
/* created: paklein (11/20/1996) */
#include "ContinuumMaterialT.h"
#include "MaterialSupportT.h"
#include "ArrayT.h"
#include "StringT.h"

namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<ContinuumMaterialT>::fByteCopy = false;
DEFINE_TEMPLATE_STATIC const bool ArrayT<ContinuumMaterialT*>::fByteCopy = true;
} /* namespace Tahoe */

using namespace Tahoe;

/* constructor */
ContinuumMaterialT::ContinuumMaterialT(void):
	ParameterInterfaceT("continuum_material"),
	fMaterialSupport(NULL),
	fNumDOF(0),
	fNumSD(0),
	fNumIP(0)
{

}

/* set the material support or pass NULL to clear */
void ContinuumMaterialT::SetMaterialSupport(const MaterialSupportT* support)
{
	fMaterialSupport = support;
	if (fMaterialSupport) {
		fNumDOF = fMaterialSupport->NumDOF();
		fNumSD = fMaterialSupport->NumSD();
		fNumIP = fMaterialSupport->NumIP();
	}
	else {
		fNumDOF = 0;
		fNumSD = 0;
		fNumIP = 0;
	}
}

/* destructor */
ContinuumMaterialT::~ContinuumMaterialT(void) { }

/* number of element nodes in the host element group */
int ContinuumMaterialT::NumElementNodes() const {
	ElementCardT& card = ElementCard(0);
	return card.NodesU().Length();
}

/* element card data */
int ContinuumMaterialT::NumElements(void) const
{
	return MaterialSupport().NumElements();
}

int ContinuumMaterialT::CurrElementNumber(void) const
{
	return MaterialSupport().CurrElementNumber();
}

ElementCardT& ContinuumMaterialT::ElementCard(int card) const
{
	ElementCardT* the_card = MaterialSupport().ElementCard(card);
	if (!the_card) ExceptionT::GeneralFail("ContinuumMaterialT::ElementCard");
	return *the_card;
}

ElementCardT& ContinuumMaterialT::CurrentElement(void) const
{
	ElementCardT* the_card = MaterialSupport().CurrentElement();
	if (!the_card) ExceptionT::GeneralFail("ContinuumMaterialT::CurrentElement");
	return *the_card;
}

/* storage initialization */
bool ContinuumMaterialT::NeedsPointInitialization(void) const { return false; }
void ContinuumMaterialT::PointInitialize(void) { /* nothing to do */ }

/* form of tangent matrix */
GlobalT::SystemTypeT ContinuumMaterialT::TangentType(void) const
{
	/* symmetric by default */
	return GlobalT::kSymmetric;
}

/* apply pre-conditions at the current time step */
void ContinuumMaterialT::InitStep(void) { }

/* finalize the current time step */
void ContinuumMaterialT::CloseStep(void) { }

/* update/reset internal variables */
void ContinuumMaterialT::UpdateHistory(void) { }
void ContinuumMaterialT::ResetHistory(void) { }

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int ContinuumMaterialT::NumOutputVariables(void) const { return 0; }
void ContinuumMaterialT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Free();	
}
void ContinuumMaterialT::ComputeOutput(dArrayT& output)
{
#pragma unused(output)
}

/* returns true if two materials have compatible output variables */
bool ContinuumMaterialT::CompatibleOutput(const ContinuumMaterialT& m1, 
	const ContinuumMaterialT& m2)
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
