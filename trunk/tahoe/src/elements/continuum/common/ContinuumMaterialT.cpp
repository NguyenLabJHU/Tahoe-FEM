/* $Id: ContinuumMaterialT.cpp,v 1.7 2002-11-14 17:06:39 paklein Exp $ */
/* created: paklein (11/20/1996) */
#include "ContinuumMaterialT.h"
#include "MaterialSupportT.h"
#include "ArrayT.h"
#include "StringT.h"

using namespace Tahoe;

/* constructor */
ContinuumMaterialT::ContinuumMaterialT(const MaterialSupportT& support):
	fMaterialSupport(support),
	fNumDOF(support.NumDOF()),
	fNumSD(support.NumSD()),
	fNumIP(support.NumIP())
{

}

/* destructor */
ContinuumMaterialT::~ContinuumMaterialT(void) { }

/* element card data */
int ContinuumMaterialT::NumElements(void) const
{
	return fMaterialSupport.NumElements();
}

int ContinuumMaterialT::CurrElementNumber(void) const
{
	return fMaterialSupport.CurrElementNumber();
}

ElementCardT& ContinuumMaterialT::ElementCard(int card) const
{
	ElementCardT* the_card = fMaterialSupport.ElementCard(card);
	if (!the_card) {
		cout << "\n ContinuumMaterialT::ElementCard: not available" << endl;
		throw ExceptionT::kGeneralFail;
	}
	return *the_card;
}

ElementCardT& ContinuumMaterialT::CurrentElement(void) const
{
	ElementCardT* the_card = fMaterialSupport.CurrentElement();
	if (!the_card) {
		cout << "\n ContinuumMaterialT::CurrentElement: not available" << endl;
		throw ExceptionT::kGeneralFail;
	}
	return *the_card;
}

/* initialization */
void ContinuumMaterialT::Initialize(void)
{
/* do nothing */
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

/* print parameters */
void ContinuumMaterialT::Print(ostream& out) const
{
	PrintName(out);
}

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

/***********************************************************************
* Protected
***********************************************************************/

void ContinuumMaterialT::PrintName(ostream& out) const
{
	out << " Material name:\n";
}
